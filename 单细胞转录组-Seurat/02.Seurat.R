message(crayon::bgCyan("=============== Package Version ==============="))
message("Seurat: ", packageVersion("Seurat"))
message("SeuratObject: ", packageVersion("SeuratObject"))

# =============== 参数===============
hue <- hue5
min_cells <- 3
min_features <- 200

# =============== 导入数据 ===============
WORKDIR <- "/home/alex/Seurat"
sc_dir <- "scdata"
source("src/config.R")

sc_files <- list.files(path = sc_dir,, pattern = "h5$",recursive=T)
if (length(sc_files) == 0) {
    stop("No h5 files found in the specified directory.")
}
sc_list <- purrr::map(sc_files, ~ Seurat::Read10X_h5(file.path(sc_dir, .x)))
names(sc_list) <- stringr::str_remove(basename(sc_files), "\\.h5")
sample_names <- names(sc_list)

srt <- purrr::map2(sc_list, sample_names, create_sc) %>%
    purrr::reduce(merge)
print(dim(srt))
# [1] 28438 70818
head(srt@meta.data )

meta <- data.table::fread("meta.xls", data.table = FALSE)
cells <- colnames(srt)
srt@meta.data %<>% dplyr::full_join(meta,by="orig.ident")
rownames(srt@meta.data) <- cells
srt %<>% subset(cells=rownames(srt@meta.data)[!is.na(srt$Group)])
cell_num_raw <- table(srt$Sample)
table(srt$Sample) %>% as.data.frame() %>% fwrite(file.path(outdir, "01_QC", "cell_num_raw.xls"), sep = "\t", col.names = F)

srt$Group %>% table()


s1 <- table(srt$Sample)
srt %<>% SeuratObject::JoinLayers(layers='counts')
srt %<>% scfilter(min.cells=min_cells,min.features=min_features)
# 去除双细胞
# 重新进行标准化
srt %<>% Seurat::NormalizeData()

srt %<>% Seurat::as.SingleCellExperiment() %>%
    scDblFinder::scDblFinder(samples = "Sample") %>%
    SeuratObject::as.Seurat() %>%
    subset(subset = scDblFinder.class == "singlet")
s2 <- table(srt$Sample)

filter_summary(before=s1,after=s2, outdir=file.path(outdir,'01_QC'), name='DoubletSummary')

# =============== 质控 ===============
srt %<>% scqc(species = "human",w=40,h=16, groupby = "Sample", cell_cycle = "hsa", hue = hue, outdir = file.path(outdir, "01_QC"), nFeature = c(400, 7000), nCount = 25000, mito = 20)
s3 <- table(srt$Sample)
filter_summary(before=s2,after=s3, outdir=file.path(outdir,'01_QC'), name='QualitySummary')
saveRDS(srt,file="srt.rds")

# =============== 整合&分群 ===============
srt %<>% scone(outdir = file.path(outdir, "02_Cluster"),h=12,verbose = TRUE, hue=hue, batch_key = "Sample",nfeatures =2000,resolution = 0.8, ndims = 20)
saveRDS(srt,file=file.path(outdir, "02_Cluster","Seurat,rds"))

avg_per_cluster <- Seurat::AverageExpression(srt, assays = "RNA", group.by = "seurat_clusters", slot = "data")$RNA
colnames(avg_per_cluster) <- gsub("g\\.", "", colnames(avg_per_cluster), perl = T)
data.table::as.data.table(avg_per_cluster,keep.rownames='Gene') %>% 
    data.table::fwrite(file.path(outdir, "02_Cluster", "avgExpression_cluster.xls"), sep = "\t")

# =============== find marker ===============
markers <- scdea(srt, outdir = file.path(outdir, "03_Marker"), hue = hue)
markers %>% dplyr::filter(p_val<0.05, pct.1>0.1, avg_log2FC>0.585) %>% fwrite(file.path(outdir, "03_Marker","marker.csv"))
scplot(srt, bunch = "seurat_clusters", hue = hue, outdir = file.path(outdir, "03_Marker"), marker = markers)

# =============== to loupe ===============
adata <- srt
barcodes <- fread('src/translation_3M-february-2018.txt.gz',header=F) %>% distinct()
adata %<>% SeuratObject::RenameCells(new.names=sample(barcodes$V1,size=ncol(adata)))
library(loupeR)
loupeR::create_loupe(
    count_mat = SeuratObject::GetAssayData(adata, assay = "RNA", layer = "counts"),
    clusters = list(Sample = as.factor(adata$Sample %>% as.character()), Cluster = as.factor(adata$seurat_clusters %>% as.character()),Group=as.factor(adata$Group %>% as.character())),
    projections = list(umap = SeuratObject::Embeddings(adata, "umap")),
    output_name = "Seurat",
    output_dir = file.path(outdir, "03_Marker"),
    force = TRUE,
    executable_path='src/louper-linux-x64'
)

