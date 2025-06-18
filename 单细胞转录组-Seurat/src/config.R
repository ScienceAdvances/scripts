Sys.setenv(LANGUAGE = "en")
options(
    bitmapType = "cairo",
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE,
    warnPartialMatchDollar = TRUE,
    timeout = 999999,
    R_MAX_VSIZE = 9e+18,
    future.globals.maxSize = 9e+16,
    Seurat.object.assay.version = "v5",
    repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN"),
    BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor"
)

# pkgs <- c("junjunlab/GseaVis",'immunogenomics/presto','genecell/COSGR','Vik2025/using',"ludvigla/semla","anndata","harmony","scCustomize","DropletUtils","SoupX","fields",'scDblFinder','AUCell','GO.db') %>% 
# pak::pkg_install('monocle')
# BiocManager::install('monocle')
# pak::pkg_install("dittoSeq")
# loupeR::setup()
# remotes::install_git('https://kkgithub.com/immunogenomics/presto')
# remotes::install_git('https://kkgithub.com/10XGenomics/loupeR')
# install.packages('/home/data/wd/Vik/CODEHUB/APP/loupeR_Linux.tar.gz',repo=NULL)
# remotes::install_local('BPCells/r')
using::using(
    jsonlite, magick, semla, Seurat, data.table, tidyverse, magrittr, using,
    conflicted, crayon, glmGamPoi, patchwork, ggpubr, anndata, reticulate,
    COSG, harmony, ggsci, future, clustree, cowplot, AUCell, pheatmap,scDblFinder,
    RColorBrewer, paletteer, scales, optparse, clusterProfiler, reshape2, scCustomize,dittoSeq
)
future::plan("future::multisession", workers = 4)
hue1 <- c("#6699A1", "#A73D7C", "#1F3683", "#A6A6A6", "#F18D8D", "#E5BA88", "#86CEEB", "#59417F")
hue2 <- c("#6699A1", "white", "#A73D7C")
hue3 <- using::hue("NPG")
hue5 <- c(
    "#a2d2e7", "#67a8cd", "#ffc17f", "#cf9f88", "#6fb3a8", "#b3e19b", "#50aa4b",
    "#ff9d9f", "#f36569", "#3581b7", "#cdb6da", "#704ba3", "#9a7fbd", "#dba9a8", "#e43030", "#e99b78", "#ff8831",
    "#e75b58", "#e49ac3", "#ab3181", "#20452e", "#bb9369", "#8c529a", "#e4d1dc", "#52a65e", "#f0b971", "#f2b09e",
    "#d4e6a1", "#55c0f2", "#496d87", hue3
) %>% unique()
dot_hue <- RColorBrewer::brewer.pal(9, 'RdBu') %>% rev()
source("src/scone.R")
source("src/scdea.R")
source("src/scqc.R")
source("src/scplot.R")
source("src/marker_dotplot.R")
source("src/utils.R")
source("src/pseudobulk_deseq.R")
source("src/pseudobulk_edger.R")
source("src/enrich.R")


outdir <- file.path(WORKDIR, "Result")
using::mkdir(outdir)
using::mkdir(file.path(outdir, "01_QC"))
using::mkdir(file.path(outdir, "02_Cluster"))
using::mkdir(file.path(outdir, "03_Marker"))
using::mkdir(file.path(outdir, "04_Annotation"))
using::mkdir(file.path(outdir, "05_Subset"))
using::mkdir(file.path(outdir, "06_Reflection"))
using::mkdir(file.path(outdir, "07_Pseudobulk"))
using::mkdir(file.path(outdir, "08_Monocle2"))
using::mkdir(file.path(outdir, "09_CellChat"))
using::mkdir(file.path(outdir, "10_Scenic"))