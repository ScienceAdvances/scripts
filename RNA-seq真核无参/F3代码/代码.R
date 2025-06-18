# remotes::install_github("junjunlab/ClusterGVis")
source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(data.table,tidyverse,ClusterGVis,ComplexHeatmap)
hue_day <- c('#F6E0D2','#DFA398','#9C6755','#659794','#EA967C','#F5C98E','#D65B5A','#586085')

exps=fread('F2/数据/salmon.gene.TPM.not_cross_norm') %>% column_to_rownames('V1')
exps=log2(exps+1)

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)


genes <- list(g1,g2,g3) %>% purrr::reduce(intersect)
exps <- exps[genes,]
dim(exps)
meta <- data.frame(Group=rep(c('CK5-6','CK7-8','CK9-10','NPA5-6', "NPA7-8", "NPA9-10"),each=3),row.names=colnames(exps))
set.seed(5201314)
exps[1:3,1:3]
pdf('F3/Cluster.pdf',width=9,height=9)
ClusterGVis::getClusters(exp = exps)
dev.off()
cm <- clusterData(exp = exps,cluster.method = "mfuzz",cluster.num = 8,seed=5201314)
# change color
pdf('F3/line.pdf',width=20,height=60)
ClusterGVis::visCluster(object = cm, ncol=1,add.mline=T, plot.type = "line") # , median, 
dev.off()


# group info
Group = meta$Group
names(Group) <- rownames(meta)
# TOP annotations
HeatmapAnnotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = Group,
    col = list(Group = c(
        "CK5-6" = hue_day[1],
        "CK7-8" = hue_day[2],
        "CK9-10" = hue_day[3],
        "NPA5-6" = hue_day[4],
        "NPA7-8" = hue_day[5],
        "NPA9-10" = hue_day[6]
    )),
    gp = grid::gpar(col = NA)
    )

# change anno bar color
pdf('F3/A.pdf',width=10,height=10)
visCluster(object = cm,
plot.type = "both", 
column_names_rot = 45,
add.box = F, 
ctAnno.col = hue_day,
HeatmapAnnotation = HeatmapAnnotation,
sample.group = meta$Group
)
dev.off()

fwrite(cm$long.res,'F3/Cluster.csv')

# 注释
d=fread('F3/Cluster.csv')
source('Code/enrich.R')
for(i in seq(8)){
    d %>% dplyr::filter(cluster==i) %>% pull('gene') %>% unique() %>% 
    enrich(suffix=str_glue("_Cluster_{i}"),outdir='F3/富集分析')
}

# KEGG热图
kegg=map_dfr(list.files('F3/富集分析',pattern='Table_KEGG'),function(x){
    cluster=str_remove_all(x,'(Table_KEGG__)|(.xlsx)|(luster)|(_)')
    d=readxl::read_excel(file.path('F3/富集分析',x))
    d$Cluster=cluster
    return(d)
})
head(kegg,1)
k2 <- kegg %>% dplyr::filter(pvalue<0.01) %>% dplyr::distinct(Description,Cluster,.keep_all = T)
mat <- k2 %>% dplyr::select(Description,p.adjust,Cluster) %>% dplyr::arrange(Cluster) %>% 
    tidyr::pivot_wider(id_cols=Description,names_from=Cluster,values_from=p.adjust) %>% 
    column_to_rownames('Description') %>% as.matrix()
using::using(ComplexHeatmap)
pdf('F3/KEGG.pdf',width = 9,height = 12)
ComplexHeatmap::Heatmap(
    mat,
    border=T,
    cluster_rows = F,
    cluster_columns = F,
    row_names_side='left',
    column_names_side = "top",
    show_column_names = TRUE,
    show_row_names = TRUE,
    column_names_rot = 0,
    width = 2, 
    height = 2, 
    na_col = "white",
    rect_gp=grid::gpar(col="black",lwd=0.8),#色块黑色色描边
    col = colorRampPalette(c('#DF1414','white'))(100),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (is.na(mat[i, j])) {
          label <- ''
        } else {
          label <- format(mat[i, j], scientific = TRUE, digits = 3)
        }
        grid::grid.text(label, x, y, gp = gpar(fontsize = 10))
      }
    )
dev.off()

# GO
go=map_dfr(list.files('F3/富集分析',pattern='Table_GO'),function(x){
    cluster=str_remove_all(tools::file_path_sans_ext(x),'Table_GO[A-Z]{2}__Cluster_') %>% paste0('C',.)
    d=readxl::read_excel(file.path('F3/富集分析',x))
    d$Cluster=cluster
    d$Type=str_split(x,'_',simplify=T)[,2]
    return(d)
})

g2 <- go %>% dplyr::filter(p.adjust<0.05) %>% dplyr::distinct(Description,Cluster,Type,.keep_all = T)
mat <- g2 %>% dplyr::select(Description,p.adjust,Cluster,Type) %>% dplyr::arrange(Cluster) %>% 
    tidyr::pivot_wider(id_cols=Description,names_from=Cluster,values_from=p.adjust) %>% 
    column_to_rownames('Description') %>% as.matrix()
g5 <- g2  %>% dplyr::distinct(Description,Type) %>% column_to_rownames('Description')

pdf('F3/GO.pdf',width = 9,height = 22)
ComplexHeatmap::Heatmap(
    mat,
    border=T,
    cluster_rows = F,
    cluster_columns = F,
    row_names_side='left',
    column_names_side = "top",
    show_column_names = TRUE,
    show_row_names = TRUE,
    column_names_rot = 0,
    width = 2, 
    height = 2, 
    right_annotation=ComplexHeatmap::HeatmapAnnotation(
        Type=g5[rownames(mat),'Type'],which='row',
        col = list(Type=c(GOBP=hue_day[1],GOMF=hue_day[3],GOCC=hue_day[4]),
        gp = gpar(col = "red", fontsize = 12)
        )),
    na_col = "white",
    rect_gp=gpar(col="black",lwd=0.8),#色块黑色色描边
    col = colorRampPalette(c('#DF1414','white'))(100),
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (is.na(mat[i, j])) {
          label <- ''
        } else {
          label <- format(mat[i, j], scientific = TRUE, digits = 3)
        }
        grid.text(label, x, y, gp = gpar(fontsize = 10))
      }
    )
dev.off()
BiocManager::install('ComplexHeatmap')
packageVersion('ComplexHeatmap')
remove.packages('ComplexHeatmap')
install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/ComplexHeatmap_2.22.0.tar.gz',repos = NULL)