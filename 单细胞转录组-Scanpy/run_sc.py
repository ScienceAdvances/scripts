# mamba activate SC
# ipython
import argparse
parse = argparse.ArgumentParser(description="scRNA-seq")
parse.add_argument("--workdir", help="workdir", type=str, required=True)
parse.add_argument("--outdir", help="outdir", type=str, required=True)
parse.add_argument("--meta", help="meta file path", type=str, required=True)
parse.add_argument(
    "--species",
    help="species ABBR. one of hsa, mmu, dre, rno",
    default="hsa",
    required=True,
)
parse.add_argument("--mt_pattern", help="mt_pattern",default=None, type=str, required=False)
parse.add_argument("--mt_gene", help="mt_gene",default=None, type=str, required=False)
parse.add_argument("--resolution", help="resolution",default=0.8, type=float, required=False)
# arg =["--workdir=/home/zktbk0/Alex/Scanpy","--outdir=/home/zktbk0/Alex/Scanpy/Result",
#        "--meta=meta.tsv","--species=hsa","--resolution=0.8","--mt_pattern=^MT-"
#        ]
# args = parse.parse_args(arg)
args = parse.parse_args()
meta = args.meta
species = args.species
resolution = args.resolution
workdir = args.workdir
outdir = args.outdir
mt_pattern = args.mt_pattern
mt_gene = args.mt_gene

import os
from pathlib import Path
import re
import warnings
import numpy as np
from scipy import sparse
import scipy.stats as spss
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import atopos
import cellscope as cs

atopos.tl.mkdir(outdir)

# setting
warnings.filterwarnings(action="ignore")
sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white')
plt.style.use('seaborn-v0_8-white')
saveimg = atopos.tl.saveimg(formats=('pdf','png'),outdir=outdir)

sc_metas=pd.read_csv(meta, sep='\t',index_col='Sample')
adata_dict = dict()
for row in sc_metas.itertuples():
    a=sc.read_10x_h5(row.h5path)
    a.obs['Group']=row.Group
    a.obs_names_make_unique()
    a.var_names_make_unique()
    adata_dict[row.Index]=a

adata = sc.concat(adata_dict,label="Sample")
adata.write_h5ad(Path(outdir).joinpath('adata_raw.h5ad'))
# adata=sc.read_h5ad(Path(outdir).joinpath('adata_raw.h5ad'))
print(">>> before filtering ncells = ",adata.n_obs)

adata.var_names_make_unique()
adata.obs_names_make_unique()

sc.pp.scrublet(adata,batch_key="Sample",expected_doublet_rate=0.05)
sc.pl.scrublet_score_distribution(adata,show=False)
atopos.tl.mkdir(Path(outdir).joinpath('A_fastqc'))
saveimg(Path('A_fastqc').joinpath("scrublet"))

adata=adata[~adata.obs.predicted_doublet]
print(">>> after doublet removing ncells = ",adata.n_obs)
cs.pp.flag_gene_family(adata,species=species)
cs.pp.flag_gene_family(adata,gene_family_name='Mito',gene_family_pattern=mt_pattern,gene_list=mt_gene)
cs.pp.fastqc(adata,
    sample="Sample",
    percent_top=(50,),
    outdir=Path(outdir).joinpath('A_fastqc'),
    qc_vars=['Mito'],
    min_genes=500,
    min_cells=10)
adata = cs.pp.mad_filter(
    adata,
    ("log1p_total_counts", 3),
    ("log1p_n_genes_by_counts", 3),
    ("pct_counts_in_top_50_genes", 3),
)
adata = adata[adata.obs.pct_counts_Mito<20,]
print(">>> after fastqc ncells = ",adata.n_obs)

# 质控后图
cs.pp.fastqc(
    adata,
    sample="Sample",
    outdir=Path(outdir).joinpath('B_qcplot_after'),
    qc_vars=['Mito'],
)

adata = cs.pp.normalise(
    adata,
    batch_key="Sample",
    n_top_genes=3000,
    outdir=Path(outdir).joinpath('C_Normalise'),
    n_jobs=6,
    regress_out=None,
    scale=False,
)

cs.pl.plot_batch_effect(
    adata,
    cluster_key="Cluster_before_integratation",
    outdir=Path(outdir).joinpath('D_batch_effect_before_integratation'), 
    batch_key="Sample",
    neighbors_key="X_pca",
    use_rep="X_pca",
)

sc.external.pp.harmony_integrate(
    adata, key="Sample", max_iter_harmony=20, adjusted_basis="X_harmony"
)

cs.pl.plot_batch_effect(
    adata,
    cluster_key="Cluster",
    outdir=Path(outdir).joinpath('E_batch_effect_after_integratation'),
    batch_key="Sample",
    neighbors_key="X_harmony",
    use_rep="X_harmony",
)

adata.write_h5ad(Path(outdir).joinpath('adata.h5ad'), compression="lzf")

atopos.tl.mkdir(Path(outdir).joinpath("F_FindAllMarkers"))
all_markers=cs.tl.find_all_markers(adata, groupby='Cluster', use_raw = False)
all_markers.to_csv(Path(outdir).joinpath("F_FindAllMarkers",'AllMakers.xls'),sep='\t')
all_markers.loc[np.logical_and(all_markers.PTS.values>0.1 , all_markers.Pvalue.values<0.05) ,:].to_excel(Path(outdir).joinpath("F_FindAllMarkers",'Makers.xlsx'))

sc.pl.rank_genes_groups(adata,key='AllMakers',show=False)
saveimg(Path('F_FindAllMarkers').joinpath("rank_genes_groups"))

sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True,key='AllMakers',show=False)
saveimg(Path('F_FindAllMarkers').joinpath("rank_genes_groups_heatmap"))

sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=5,key='AllMakers',show=False)
saveimg(Path('F_FindAllMarkers').joinpath("rank_genes_groups_stacked_violin"))

sc.pl.rank_genes_groups_matrixplot(adata,n_genes=5,key='AllMakers',show=False)
saveimg(Path('F_FindAllMarkers').joinpath("rank_genes_groups_matrixplot"))
