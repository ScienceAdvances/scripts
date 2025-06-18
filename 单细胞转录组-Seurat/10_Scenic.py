#!/usr/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/
mamba activate scenic
mamba activate R
source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(anndataR, loomR,Seurat,magrittr)
pak::pak("scverse/anndataR", dependencies = TRUE)
remotes::install_github('mojaveazure/loomR')

srt=readRDS('结果/03.Annotation/Seurat.rds')
srt %<>% SeuratObject::UpdateSeuratObject()
srt[["RNA"]] <- as(object = srt[["RNA"]], Class = "Assay5")
srt[["HTO"]] <- as(object = srt[["HTO"]], Class = "Assay5")
adata=srt
adata$CellType %<>% as.character()
adata@meta.data[adata$CellType=='Naïve CD8+ T cells','CellType'] = 'Naive CD8+ T cells'
anndataR::write_h5ad(adata,path='结果/Scenic/adata.h5ad',compression='lzf')
loomR::create(
    filename = "结果/Scenic/X.loom",
    overwrite=TRUE,
    gene.attrs=list(var_names=rownames(adata)),
    cell.attrs=list(obs_names=colnames(adata)), 
    data=SeuratObject::GetAssayData(adata,layer='data',assay='RNA'),
    do.transpose = TRUE)
saveRDS(adata,file='结果/03.Annotation/Seurat.rds')

mamba activate scenic
n_jobs=12
databse=$HOME/Vik/S134/SCENIC
workdir=/home/data/wd/Vik/S134/结果/Scenic
tfs_fname=hs_hgnc_tfs.txt
feather=hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
annotations_fname=motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
tfs_fname=allTFs_mm.txt
feather=mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
annotations_fname=motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
expression_mtx_fname=X.loom


# STEP 1/3: Gene regulatory network inference, and generation of co-expression modules
pyscenic grn \
    --method grnboost2 \
    --num_workers $n_jobs \
    --seed 1314520 \
    --sparse \
    --cell_id_attribute obs_names \
    --gene_attribute var_names \
    --output $workdir/grn.csv \
    $workdir/$expression_mtx_fname \
    $databse/$tfs_fname
# STEP 2/3: Regulon prediction (cisTarget)
pyscenic ctx \
    --mask_dropouts \
    --num_workers $n_jobs \
    --sparse \
    --cell_id_attribute obs_names \
    --gene_attribute var_names \
    --annotations_fname $databse/$annotations_fname \
    --expression_mtx_fname $workdir/$expression_mtx_fname \
    --output $workdir/ctx.csv \
    $workdir/grn.csv \
    $databse/$feather



def filter_regulons(
    motifs,
    db_names=("mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings",),
):
    """
    从ctx.csv中筛选重要的regulons，之后再运行AUCell
    """
    motif = df_motifs.copy()
    motif.columns = motif.columns.droplevel(0)

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motif for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.

    context = [".".join(list(x)) for x in motif.Context.values]

    b1=bq.st.detects(string=context,pattern=r"{}".format("|".join(db_names)))
    b2=bq.st.detects(string=context,pattern="weight>50.0%")
    b3=bq.st.detects(string=context,pattern="activating")
    motif=motif.loc[b1 & ~b2 & b3, :]

    # We build regulons only using enriched motif with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    b1 = motif.NES.values >= 3.0
    b2 = motif.Annotation.values == "gene is directly annotated"
    b3 = motif.Annotation.str.startswith("gene is orthologous to").values
    b4 = motif.Annotation.str.endswith("which is directly annotated for motif").values
    regulons = list(
        filter(lambda r: len(r) >= 10, df2regulons(motif.loc[b1 & (b2 | b3 & b4), :]))
    )

    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))
OUTPUT_DIR='/home/data/wd/Vik/S134/结果/Scenic/PBS'
outdir='/home/data/wd/Vik/S134/结果/Scenic'

adata = sc.read_h5ad('结果/Scenic/adata.h5ad')

adata=adata[adata.obs.Group=='PBS',:]
df_motifs = load_motifs(f"{outdir}/ctx.csv")
regulons = filter_regulons(df_motifs)
adata.layers['counts'].shape
adata.X = None
adata.X = adata.layers['data']
auc_mtx = aucell(exp_mtx=adata.to_df(), signatures=regulons, seed=1314520, num_workers=12)
auc_mtx.to_csv(f"{OUTPUT_DIR}/aucell.csv")
np.sum(auc_mtx>0,axis=1)
bin_mtx, thresholds = binarize(auc_mtx,seed=1314520,num_workers=12)
bin_mtx.to_csv(f"{OUTPUT_DIR}/bin_mtx.csv")
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(f"{OUTPUT_DIR}/thresholds.csv")

auc_mtx=pd.read_csv(OUTPUT_DIR+"/aucell.csv", index_col=0)
bin_mtx = pd.read_csv(OUTPUT_DIR+"/bin_mtx.csv", index_col=0)
thresholds = pd.read_csv(OUTPUT_DIR+"/thresholds.csv", index_col=0).threshold
# 删除基因后的(+)
auc_mtx.columns = bq.st.removes(string=auc_mtx.columns, pattern=r'\(\+\)')
bin_mtx.columns = bq.st.removes(string=bin_mtx.columns, pattern=r'\(\+\)')
thresholds.index = bq.st.removes(string=thresholds.index, pattern=r'\(\+\)')

auc_sum = auc_mtx.apply(sum,axis=0).sort_values(ascending=False)
fig, axes = plt.subplots(1, 5, figsize=(8, 2), dpi=300)
for x,y in enumerate(axes):
    plot_binarization(auc_mtx, auc_sum.index[x], thresholds[auc_sum.index[x]], ax=y)
plt.tight_layout()
plt.savefig(OUTPUT_DIR+"/plot_binarization.pdf",bbox_inches='tight')

cell_type_key = "CellType"
cell_type_color_lut = dict(zip(adata.obs[cell_type_key].unique(), sc.pl.palettes.default_20))
cell_id2cell_type_lut = adata.obs[cell_type_key].to_dict()
bw_palette = sns.xkcd_palette(["white", "black"])
sns.set()
sns.set(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})

g = sns.clustermap(
    data=bin_mtx.T.loc[:,adata.obs_names], 
    col_colors=auc_mtx.loc[adata.obs_names,:].index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
    cmap=bw_palette, 
    figsize=(20,20)
    )
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)
plt.savefig(outdir+"/clustermap.pdf",bbox_inches='tight')

# logo
def fetch_logo(regulon, base_url = "http://motifcollections.aertslab.org/v10/logos/"):
    for elem in regulon.context:
        if elem.endswith('.png'):
            return '<img src="{}{}" style="max-height:124px;"></img>'.format(base_url, elem)
    return ""
df_regulons = pd.DataFrame(data=[list(map(operator.attrgetter('name'), regulons)),
    list(map(len, regulons)),
    list(map(fetch_logo, regulons))], 
    index=['name', 'count', 'logo']).T
df_regulons.to_csv(outdir+"/logo.csv",index=False)
with open(f'{outdir}/logo.html',mode='w') as f:
    f.write(IPython.display.HTML(df_regulons.to_html(escape=False)).data)

# 细胞特异 REGULATORS
from pyscenic.export import add_scenic_metadata
add_scenic_metadata(adata, auc_mtx.loc[adata.obs_names,:], regulons);
df_scores=bq.tl.select(adata.obs,columns=[cell_type_key],pattern=r'^Regulon\(')
df_results = ((df_scores.groupby(by=cell_type_key).mean() - df_scores[df_scores.columns[1:]].mean())/ df_scores[df_scores.columns[1:]].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
df_results['Regulon'] = bq.st.removes(string=df_results.regulon,pattern='(Regulon\()|(\)_[xy])')

df_heatmap = pd.pivot_table(data=df_results[df_results.Z >= 0.5].sort_values('Z', ascending=False),index=cell_type_key, columns='Regulon', values='Z')
fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
sns.heatmap(df_heatmap, ax=ax1, annot=True, fmt=".1f", linewidths=.7, cbar=False, square=True, linecolor='gray',cmap="viridis", annot_kws={"size": 8})
ax1.set_ylabel('');
plt.savefig(outdir+"/REGULATORS_hetamap.pdf",bbox_inches='tight')


from pyscenic.rss import regulon_specificity_scores
rss = regulon_specificity_scores(auc_mtx.loc[adata.obs_names,:], adata.obs[cell_type_key])
for i in adata.obs.CellType.unique():
    plot_rss(rss, cell_type=i, top_n=5)
    plt.savefig(f"{outdir}/{i}_RSS.pdf",bbox_inches='tight')

aucell_adata = sc.AnnData(X=auc_mtx.loc[adata.obs_names,:].sort_index(),dtype=np.float32)
aucell_adata.obs = adata.obs
names = list(map(operator.attrgetter('name'), filter(lambda r: r.score > 4.0, regulons)))
sc.pl.stacked_violin(aucell_adata, names, groupby=cell_type_key)
plt.savefig(f"{outdir}/stacked_violin.pdf",bbox_inches='tight')

rss_df=rss.melt(var_name="Regulon", value_name  = 'RSS',ignore_index=False)
rss_df.loc[:,"Key"] = [f"{x}-{y}" for x,y in zip(rss_df.index,rss_df.Regulon)]
df_results.loc[:,"Key"] = [f"{x}-{y}" for x,y in zip(df_results.CellType,df_results.Regulon)]
merge_data=pd.merge(df_results,rss_df.drop(columns="Regulon"),on="Key")
plot_data=merge_data.groupby(cell_type_key).apply(lambda x: x.sort_values(by="RSS",ascending=False).iloc[0:6,:])
plot_data.loc[:,"Regulon"]=bq.st.removes(string=plot_data.Regulon,pattern=r"\(\+\)")
plot_data.drop_duplicates("Regulon",inplace=True)
plot_data.loc[:,"y"] = pd.Categorical(plot_data.Regulon,categories=plot_data.Regulon,ordered=True)
from plotnine import *
p=(ggplot(plot_data)+
    geom_point(aes(cell_type_key,"y",color = "RSS", size = "Z"))+
    labs(title="",x="",y="")+
    theme_bw()+
    theme(axis_text_x = element_text(angle = 45, hjust = 1))
    
)
ggsave(p,filename=f"{outdir}/dotplot.pdf",width=4,height=5)
merge_data.to_csv(f"{outdir}/CellType_Regulon.csv",index=False)

