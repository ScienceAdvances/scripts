rm(list = ls())
source("src/config.R")

hue3 <- colorls$NPG
options(Seurat.object.assay.version = "v3")
using::using(scDblFinder,scCustomize)
# === === === === === 导入数据 ===  === === === === ===
WORK_DIR <- "/home/alex/HTO"
outdir <- '结果/SC01'
dir.create(outdir,recursive=T)
setwd(WORK_DIR)
sc_dir <- "scdata"
min_cells <- 3
min_features <- 200
sc_files <- list.files(path = sc_dir, pattern = "h5$")
sc_list <- purrr::map(sc_files, ~ Seurat::Read10X_h5(file.path(sc_dir, .x)))
names(sc_list) <- stringr::str_remove(sc_files, "\\.h5")
sample_names <- names(sc_list)

create_sc <- function(x, y) {
  colnames(x) <- paste(colnames(x), y, sep = "_")
  srt <- Seurat::CreateSeuratObject(
    counts = x,
    min.cells = min_cells,
    min.features = min_features,
    project = y
  )
  return(srt)
}

srt <- purrr::map2(sc_list, sample_names, create_sc) %>%
  purrr::reduce(merge)

mgs=Seurat::Read10X_h5('scdata/MGS/outs/filtered_feature_bc_matrix.h5')
pbs=Seurat::Read10X_h5('scdata/PBS/outs/filtered_feature_bc_matrix.h5')
mgs=CreateAssayObject(mgs[['Antibody Capture']])
colnames(mgs) %<>% paste('MGS', sep = "_")
pbs=CreateAssayObject(pbs[['Antibody Capture']])
colnames(pbs) %<>% paste('PBS', sep = "_")
hto <- merge(mgs,pbs)
common_cells <- base::intersect(colnames(hto),colnames(srt))
srt %<>% subset(cells=common_cells)
hto %<>% subset(cells=common_cells)
srt[['HTO']] <- hto

print(dim(srt))
# [1] 18412 16030

srt@meta.data$orig.ident %>% table()
#  MGS  PBS 
# 8240 7790 
DefaultAssay(srt) <- 'HTO'
srt %<>% Seurat::NormalizeData(assay = "HTO", normalization.method = "CLR")
srt %<>% Seurat::HTODemux(assay = "HTO", positive.quantile = 0.99)
table(srt$HTO_classification.global)
#  Doublet Negative  Singlet 
#     2990     1196    11844 
head(srt@meta.data,2)
# Group cells based on the max HTO signal
Idents(srt) <- "HTO_maxID"
RidgePlot(srt, assay = "HTO", features = rownames(srt[["HTO"]]), ncol = 4) %>% 
  gs(outdir=outdir,name="_RidgePlot_HTO",w = 16,h = 6)
FeatureScatter(srt, feature1 = "hto_tag1", feature2 = "hto_tag2",col=hue) %>% 
  gs(outdir=outdir,name="_Scatter_HTO",w = 8,h = 8)
Idents(srt) <- "HTO_classification.global"
VlnPlot(srt, features = "nCount_RNA", pt.size = 0.1, log = TRUE) %>% 
  gs(outdir=outdir,name="_nCount_RNA_HTO",w = 8,h = 8)

# First, we will remove negative cells from the object
srt.subset <- subset(srt, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(srt.subset) <- "HTO"
srt.subset <- ScaleData(srt.subset, features = rownames(srt.subset),verbose = FALSE)
srt.subset <- RunPCA(srt.subset, features = rownames(srt.subset), approx = FALSE)
srt.subset <- RunTSNE(srt.subset, dims = 1:4, perplexity = 100)
DimPlot(srt.subset) %>% 
  gs(outdir=outdir,name="_Doublets_HTO",w = 6,h = 6)
HTOHeatmap(srt.subset, assay = "HTO", ncells = 5000) %>% 
  gs(outdir=outdir,name="_Heatmap_HTO",w = 6,h = 6)

# Extract the singlets
srt.singlet <- subset(srt.subset, idents = "Singlet")

# Select the top 1000 most variable features
srt.singlet <- FindVariableFeatures(srt.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
srt.singlet <- ScaleData(srt.singlet, features = rownames(srt.singlet))

# Run PCA
srt.singlet <- RunPCA(srt.singlet, features = rownames(srt.singlet))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
srt.singlet %<>% FindNeighbors(reduction = "pca", dims = 1:3) %>% 
  FindClusters(resolution = 0.6, verbose = FALSE) %>%
  RunTSNE(reduction = "pca", dims = 1:3) %>% 
  RunUMAP(reduction = "pca", dims = 1:3)

# Projecting singlet identities on TSNE visualization

scCustomize::DimPlot_scCustom(srt.singlet,reduction = 'umap',group.by = 'HTO_classification',colors_use=hue,figure_plot=TRUE) %>% 
gs(outdir=outdir,name='_UMAP_HTO_classification',w=8,h=8)

srt <- srt.singlet
DefaultAssay(srt) <- "RNA"
srt$Sample <- srt$orig.ident
# === === === === === 质控 ===  === === === === ===
source('/home/data/wd/Vik/CODEHUB/SingleCell/scrna_qc.R')
srt %<>% scrna_qc(species='mouse', groupby='Sample',cell_cycle='mmu', hue=hue, outdir=outdir, nFeature=c(400,7000),  nCount=25000,  mito = 20)
cell_num_raw=table(srt$Sample)
##QC ,filter low Q cell
##output experiment raw information
out_gene_umi=as.matrix(srt@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file=paste(outdir,"Basicinfo_nGene_nUMI_mito.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

cell_num_filter=table(srt$Sample)
filter_summary=cbind(cell_num_raw,cell_num_filter)
filter_summary=data.frame("cell_num_raw"=filter_summary[,1],"cell_num_filter"=filter_summary[,2])
filter_summary[nrow(filter_summary)+1,]=apply(filter_summary,2,sum)
rownames(filter_summary)[nrow(filter_summary)]="total"
filter_summary$percent=apply(filter_summary,1,function(x){round(x[2]*100/x[1],2)})

write.table(filter_summary,paste(seurat_qc_dir,"cells_filter_stat.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#for more than one sample
#split data by sample
#names(ifnb.list)=sample_list
cluster_dir <- outdir
opt <- list(
    gather="harmony",
    type='hsa'
)
sample_list <- c('MGS','PBS')

srt %<>% blend(outdir, batch_key="Sample", resolution=0.8, ndim=20, method='harmony')

saveRDS(srt,file=file.path(outdir,"seuart.rds"))

avg_per_cluster=AverageExpression(srt,"RNA")$RNA
colnames(avg_per_cluster)=gsub("RNA.","",colnames(avg_per_cluster),perl=T)
write.table(avg_per_cluster,paste(cluster_dir,"avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#find marker
markers <- findmarker(srt,outdir)
srt$Group <- srt$Sample
Seurat.Plot(srt, cluster="seurat_clusters", hue = hue3, outdir = outdir, markers = markers)


# to loupe
# install.packages('/home/data/wd/Vik/CODEHUB/APP/loupeR_Linux.tar.gz',repo=NULL)
adata=srt
# barcodes <- fread('/home/data/wd/APP/cellranger-8.0.1/lib/python/cellranger/barcodes/3M-3pgex-may-2023.txt.gz',header=F) %>% 
#     distinct()
# adata %<>% SeuratObject::RenameCells(new.names=sample(barcodes$V1,size=ncol(adata)))
loupeR::create_loupe(
  count_mat=SeuratObject::GetAssayData(adata,assay='RNA',slot='counts'),
  clusters=list(Sample=as.factor(adata$Sample),Cluster=as.factor(adata$seurat_clusters)),
  projections=list(umap=SeuratObject::Embeddings(adata,'umap')),
  output_name='S134',
  output_dir=outdir,
  force = TRUE
  )

srt <- readRDS('seuart.rds')
outdir <- 'Result/Cluster'
i=0.2
srt@meta.data %>% colnames()
srt@meta.data %<>% dplyr::select(-SCT_snn_res.0.2,-Cluster_0.2,-SCT_snn_res.0.8)
a1=srt@meta.data$seurat_clusters
a2=srt@meta.data$seurat_clusters
a1
a2
colnames()
for (i in c(0.2,0.6,0.8,1.0,1.2)) {
    cluster_name <- base::paste("Cluster", i, sep = "_")
    outdir_cluster <- file.path(outdir, cluster_name)
    using::mkdir(outdir_cluster)
    srt %<>% Seurat::FindClusters(resolution = i,graph.name='SCT_snn', verbose = TRUE)
    colnames(srt@meta.data)[colnames(srt@meta.data)=='seurat_clusters'] = cluster_name
    # Seurat::Idents(srt) <- cluster_name
    # scCustomize::DimPlot_scCustom(srt,reduction = 'umap',group.by = cluster_name,colors_use=hue3,figure_plot=TRUE) %>% 
    # gs(outdir=outdir_cluster,name=paste0('UMAP_',cluster_name),w=9,h=8)
    # scCustomize::FeaturePlot_scCustom(srt,reduction = 'umap',features = 'percent_mito',colors_use=viridis_plasma_dark_high,figure_plot=TRUE) %>% 
    # gs(outdir=outdir_cluster,name=paste0('UMAP_Mito_',cluster_name),w=9,h=8)
    # scCustomize::FeaturePlot_scCustom(srt,reduction = 'umap',features = 'percent_ribo',colors_use=viridis_plasma_dark_high,figure_plot=TRUE) %>% 
    # gs(outdir=outdir_cluster,name=paste0('UMAP_Ribo_',cluster_name),w=9,h=8)
    # markers <- findmarker(srt,outdir_cluster)
    # markers %>% dplyr::group_by(Cluster) %>% dplyr::top_n(20, wt=avg_log2FC) %>% dplyr::ungroup() %>% 
    #   dplyr::arrange(Cluster) %>% 
    #   fwrite(file=file.path(outdir_cluster,'Top20Marker.xls'),sep='\t')
    # srt$Group <- srt$Sample
    # Seurat.Plot(srt, bunch=cluster_name, hue = hue3, outdir = outdir_cluster, marker = marker)
}

clustree::clustree(srt@meta.data, prefix = "Cluster_") %>% using::gs(outdir=outdir,name='clustree',w=12,h=12)
source('src/seurat.plot.R')