# BiocManager::install("ZJU-UoE-CCW-LAB/scCDC")
start <- Sys.time()
message(crayon::bgBlue("Seurat start time: "), format(start, " %F %T"))
# 描述参数的解析方式
option_list <- list(
    optparse::make_option(c("--h5path"), type = "character", action = "store", default = FALSE, help = "h5path"),
    optparse::make_option(c("--resolution"), type = "double", action = "store", default = 0.8, help = "resolution"),
    optparse::make_option(c("--outdir"), type = "character", action = "store", default = FALSE, help = "outdir"),
    optparse::make_option(c("--sample_name"), type = "character", default = FALSE, action = "store", help = "sample_name"),
    optparse::make_option(c("--verbose"), type = "logical", default = TRUE, action = "store", help = "verbose")
)
source("/home/zktbk0/Alex/CODEHUB/config.R")
# ======================= 解析参数  =============================
opt <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list, add_help_option = TRUE,
    usage = "Usage: %prog [options] \nDescription: Remove ambient RNA!"
))
# opt <- list(
#     h5path="/data/Alex/X201/SRR11038990/outs/filtered_feature_bc_matrix.h5",
#     resolution=0.3,
#     outdir='/home/data/wd/Vik/J286-2/Result/05_MPs/Cluster',
#     sample_name='SRR11038990',
#     verbose=TRUE
# )
using::using(magrittr,scCDC,Seurat,SeuratObject,DropletUtils)
resolution <- opt$resolution
h5path <- opt$h5path
sample_name <- opt$sample_name
verbose=TRUE
outdir <- file.path(opt$outdir, sample_name)
using::mkdir(outdir)

srt = Seurat::Read10X_h5(h5path) %>% 
    Seurat::CreateSeuratObject(min.cells=3,min.features=200) %>% 
    Seurat::NormalizeData(normalization.method = "LogNormalize",scale.factor = 10000,verbose=verbose) %>% 
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000,verbose=verbose) %>% 
    Seurat::ScaleData(features = rownames(.),verbose=verbose) %>% 
    Seurat::RunPCA(features = VariableFeatures(object = .),verbose=verbose) %>% 
    Seurat::RunUMAP(dims=1:20,verbose=verbose) %>% 
    Seurat::FindNeighbors(dims = 1:20,verbose=verbose) %>%
    Seurat::FindClusters(resolution=0.8,verbose=verbose)

GCGs <- scCDC::ContaminationDetection(srt)

GCGs <- scCDC::ContaminationDetection(
    srt,
    restriction_factor = 0.5, 
    sample_name = sample_name,
    out_path.plot = paste0(outdir, "/"),
    out_path.table = paste0(outdir, "/")
)

mislet_cont_ratio <- scCDC::ContaminationQuantification(srt,rownames(GCGs))
# The maximum contamination ratio is:0.01677, which indicates high contamination in this dataset. scCDC is highly recommended to be applied.
seuratobj_corrected <- scCDC::ContaminationCorrection(srt, rownames(GCGs))

DropletUtils::write10xCounts(
    path=file.path(outdir,glue::glue("{sample_name}.h5")),
    x=SeuratObject::GetAssayData(seuratobj_corrected, layer="counts",assay="Corrected"),
    barcodes = colnames(seuratobj_corrected),
    gene.symbol = rownames(seuratobj_corrected),
    gene.id = rownames(seuratobj_corrected),
    type = "HDF5",
    overwrite = TRUE,
    version = '3'
)
# 89:0.02903
# 90:0.01677
# 91:0.03411
# 92:0.015
# 93:0.00755
# 94:0.01749
# 95:0.05101

