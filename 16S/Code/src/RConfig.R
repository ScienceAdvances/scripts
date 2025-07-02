library(using)
using::using(
    jsonlite, plyr, data.table, tidyverse, magrittr,
    conflicted, crayon, patchwork, ggpubr, reticulate,
    ggsci, future, cowplot, pheatmap,
    RColorBrewer, paletteer, optparse, clusterProfiler, reshape2,
    mecodev, phyloseq, microeco, ape, vegan, MicrobiotaProcess
)
source('/home/data/wd/Vik/X155/Code/src/functions.R')
hue4 <- c("#E54B34", "#4CBAD4", "#009F86", "#3B5387", "#F29A7F", "#8491B3", "#91D1C1", "#DC0000", "#7E6047", "#CCCCCC", "#BC8B83", "#33ADAD", "#347988", "#9F7685", "#C1969A", "#8BB0BB", "#CE8662", "#B04929", "#A59487", "#E3907E", "#D46F5B", "#41B4C1", "#278C87", "#726486", "#DA988C", "#88A0B7", "#B9AC91", "#C63517", "#927A66", "#DBAEA4", "#97A4AB", "#21A69A", "#3A6688", "#C98882", "#A593A7", "#8EC0BE", "#D85935", "#985738", "#B9AFA9", "#E67059", "#E5BFB9", "#B2CED4", "#779F99", "#747A87", "#F2DCD5", "#A7ABB3", "#C1D1CD", "#DCA5A5", "#7E7770", "#CCCCCC")
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
setheme <- function(xtext_rotation = NULL) {
    if(is.null(xtext_rotation)){
        hjust = NULL
    }else {
        hjust = 1
    }
    theme1 <- ggplot2::theme(
        axis.title = ggplot2::element_text(
            face = "bold",
            size = "18",
            color = "black"
        ),
        axis.text.x = ggplot2::element_text(
            color = "black",
            face = "bold",
            size = 10,
            angle = xtext_rotation,
            hjust = hjust
        ),
        axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black"),
        legend.text = ggplot2::element_text(face = "bold", color = "black", size = 10),
        legend.title = ggplot2::element_text(face = "bold", color = "black", size = 10),
        legend.position = "right",
        strip.text.x = ggplot2::element_text(face = "bold", size = 10),
        strip.text.y = ggplot2::element_text(face = "bold", size = 10),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = ggplot2::element_text(face = "bold", color = "black", lineheight = .8, hjust = 0.5, size = 5)
    )
    return(theme1)
}
