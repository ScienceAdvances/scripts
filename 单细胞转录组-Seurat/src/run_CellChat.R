source('src/config.R')
option_list <- list(
    make_option("--workdir", type = "character", action = "store", default = "./", help = "工作目录路径"),
    make_option("--anndata", type = "character", action = "store", default = NULL, help = "文件标识"),
    make_option("--seurat", type = "character", action = "store", default = NULL, help = "文件标识"),
    make_option("--filename", type = "character", action = "store", default = NULL, help = "文件标识"),
    make_option("--species", type = "character", action = "store", default = "hsa", help = "物种"),
    make_option("--celltype", type = "character", action = "store", default = 'CellType', help = "celltype"),
    make_option("--used_celltype", type = "character", action = "store", default = NULL, help = "used_celltype"),
    make_option("--cluster", type = "character", action = "store", default = 'Group', help = "key"),
    make_option("--subset", type = "character", action = "store", default = NULL, help = "subset"),
    make_option("--conda", type = "character", action = "store", default = '/home/data/wd/miniforge3/envs/sc', help = "conda")
    # make_option("--n_jobs", type = "integer", action = "store", default = 16, help = "线程")
)
opt <- parse_args(OptionParser(
    option_list = option_list, add_help_option = TRUE,
    usage = "Usage: %prog [options] \nDescription: CellChat pre!"
))
using::using(CellChat)
reticulate::use_condaenv(condaenv = opt$conda, required = TRUE)
# opt <- list(
#     anndata="/home/data/wd/Will/S128/adata.h5ad",
#     filename='SLE',
#     cluster='Group2',
#     celltype='CellType',
#     subset='SLE',
#     # seurat='/home/data/wd/victor/JL139/结果/reflect/reflect.seurat',
#     species='hsa',
#     workdir="/home/data/wd/Will/S128/结果/SC05",
#     conda='/home/data/wd/miniforge3/envs/sc',
#     # used_celltype='/home/data/wd/victor/JL139/结果/CellChat/used_celltype.txt'，
#     n_jobs=16
#     )
dir.create(opt$workdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$workdir)

# ==============================================================================
# 1.创建CellChat对象
# ==============================================================================
if(!is.null(opt$anndata)){
    ad <- reticulate::import("anndata")
    adata <- ad$read_h5ad(opt$anndata)
    if(!is.null(subset)){
        adata=adata[adata$obs[opt$cluster]==opt$subset,]$copy()
    }
    if(adata$n_obs > 2e+5){
        meta=adata$obs
        meta %<>% dplyr::slice_sample(prop=0.3,by=opt$celltype)
        adata=adata[rownames(meta),]$copy()
    }
    if(!is.null(opt$used_celltype)){
        used_celltype <- strsplit(opt$used_celltype,',') %>% unlist()
        meta=adata$obs
        meta %<>% dplyr::dplyr(meta[celltype] %in% used_celltype)
        adata=adata[rownames(meta),]$copy()
    }
    cc <- CellChat::createCellChat(adata$T$X, meta=adata$obs, group.by = opt$celltype)
}

if(!is.null(opt$seurat)){
    srt <- base::readRDS(opt$seurat)
    if(!is.null(opt$subset)){
        subset <- strsplit(opt$subset,',') %>% unlist()
        SeuratObject::Idents(srt) <- opt$cluster
        srt %<>% subset(idents=opt$subset)
    }
    
    SeuratObject::Idents(srt) <- opt$celltype
    celltypes <- SeuratObject::Idents(srt) %>% levels()
    if(!is.null(opt$used_celltype)){
        used_celltype <- strsplit(opt$used_celltype, ',') %>% unlist()
        SeuratObject::Idents(srt) <- opt$celltype
        srt %<>% subset(idents=used_celltype)
        SeuratObject::Idents(srt) <- "CellType"
    }
    if(ncol(srt) > 2e+5){
        meta=srt$meta.data
        meta %<>% dplyr::slice_sample(prop=0.3,by=opt$celltype)
        srt %<>% subset(cells=rownames(meta))
    }
    cc <- CellChat::createCellChat(srt, group.by = opt$celltype, assay = "RNA")
}
db_used <- switch(opt$species,
    hsa = CellChat::CellChatDB.human,
    mmu = CellChat::CellChatDB.mouse,
    dre = CellChat::CellChatDB.zebrafish
)
db_used$interaction <- db_used$interaction %>% dplyr::distinct(interaction_name_2,.keep_all = TRUE)
cc@DB <- db_used

# ==============================================================================
# 2.预处理表达数据 + 细胞通讯网络
# ==============================================================================
future::plan("multisession", workers = 1)
cc %<>% CellChat::subsetData()
cc %<>% CellChat::identifyOverExpressedGenes()
cc %<>% CellChat::identifyOverExpressedInteractions()
cc %<>% CellChat:::smoothData(adj=(switch(opt$species,
    hsa = CellChat::PPI.human,
    mmu = CellChat::PPI.mouse,
    dre = CellChat::PPI.zebrafish
)))
cc %<>% CellChat::computeCommunProb()
cc %<>% CellChat::filterCommunication(min.cells = 10)


cc %<>% CellChat::computeCommunProbPathway()
cc %<>% CellChat::aggregateNet()
cc %<>% CellChat::netAnalysis_computeCentrality(slot.name = "netP")
cc %<>% CellChat::computeNetSimilarity(type = "functional")
cc %<>% CellChat::netEmbedding(type = "functional", umap.method='umap-learn')

cc %<>% CellChat::netClustering(fig.id = opt$filename, type = "functional", do.parallel = F, nCores=1)
cc %<>% CellChat::computeNetSimilarity(type = "structural")
cc %<>% CellChat::netEmbedding(type = "structural",umap.method='umap-learn')
cc %<>% CellChat::netClustering(fig.id = opt$filename, type = "structural", do.parallel = F, nCores=1)
saveRDS(cc, file = str_glue("CellChat_{opt$filename}.rds"))
