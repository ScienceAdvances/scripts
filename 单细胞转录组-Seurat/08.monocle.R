#!/usr/bin/Rscript
using::using(crayon,monocle,anndata,patchwork,SingleCellExperiment,magrittr,stringr,Seurat,optparse,DropletUtils)
message(crayon::bgCyan("=============== Package Version ==============="))
message("Seurat: ", packageVersion("Seurat"))
message("SeuratObject: ", packageVersion("SeuratObject"))
message("monocle: ", packageVersion("monocle"))
# BiocManager::install('monocle')
start <- Sys.time()
message(crayon::bgBlue("monocle2 start time: "), format(start, " %F %T"))

option_list <- list(
    make_option(c("--seurat"), type = "character", action = "store", default = NULL, help = "seurat"),
    make_option(c("--anndata"), type = "character", action = "store", default = NULL, help = "resolution"),
    make_option(c("--outdir"), type = "character", action = "store", default = NULL, help = "outdir"),
    make_option(c("--filename"), type = "character", action = "store", default = NULL, help = "filename"),
    make_option(c("--cluster"), type = "character", action = "store", default = NULL, help = "cluster"),
    make_option(c("--celltype"), type = "character", action = "store", default = NULL, help = "celltype"),
    make_option(c("--subset"), type = "character", action = "store", default = NULL, help = "subset"),
    make_option(c("--conda"), type = "character", action = "store", default = NULL, help = "conda"),
    make_option(c("--used_celltype"), type = "character", action = "store", default = NULL, help = "used_celltype"),
    make_option(c("--n_jobs"), type = "integer", action = "store", default = 2, help = "n_jobs")
)

# ======================= 解析参数  =============================
opt <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list, add_help_option = TRUE,
    usage = "Usage: %prog [options] \nDescription: Remove ambient RNA!"
))
# opt <- list(
#     seurat='/home/data/wd/Vik/S147/Result/06_Macrophages/Annotation/seurat.rds',
#     # anndata="/home/data/wd/Vik/S128/结果/SC06/adata.h5ad",
#     filename='Normal', # Tumor
#     cluster='Group',
#     celltype='CellType',
#     subset='Normal',
#     # seurat='/home/data/wd/victor/JL139/结果/reflect/reflect.seurat',
#     species='hsa',
#     outdir="/home/data/wd/Vik/S147/Result/07.monocle2/Normal",
#     conda='/home/data/wd/miniforge3/envs/sc',
#     # used_celltype='/home/data/wd/victor/JL139/结果/CellChat/used_celltype.txt'，
#     n_jobs=2
#     )
source('src/config.R')
reticulate::use_condaenv(opt$conda)
dir.create(opt$outdir,recursive = T)
setwd(opt$outdir)
# cds <- readRDS('monocle2.rds')
outdir=opt$outdir

# ==============================================================================
# 1.创建CellChat对象
# ==============================================================================
if(!is.null(opt$anndata)){
    ad <- reticulate::import("anndata")
    adata <- ad$read_h5ad(opt$anndata)
    if(!is.null(subset)){
        adata=adata[adata$obs[opt$cluster]==opt$subset,]$copy()
    }
    if(adata$n_obs > 2e+4){
        meta=adata$obs
        meta %<>% dplyr::slice_sample(prop=0.3,by=c(opt$celltype,opt$cluster))
        adata=adata[rownames(meta),]$copy()
    }
    if(!is.null(opt$used_celltype)){
        used_celltype <- split(opt$used_celltype,',') %>% unlist()
        meta=adata$obs
        meta %<>% dplyr::dplyr(meta[celltype] %in% used_celltype)
        adata=adata[rownames(meta),]$copy()
    }
    fdata <- adata$var
    fdata$gene_short_name <- rownames(fdata)
    cds <- monocle::newCellDataSet(
        cellData=Matrix::Matrix(adata$T$X),
        phenoData = new("AnnotatedDataFrame", data = adata$obs),
        featureData = new("AnnotatedDataFrame", data = fdata),
        expressionFamily = VGAM::negbinomial.size()
    )
    hvg <- rownames(fdata)[fdata$highly_variable]
}

if(!is.null(opt$seurat)){
    srt <- base::readRDS(opt$seurat)
    if(!is.null(opt$subset)){
        SeuratObject::Idents(srt) <- opt$cluster
        srt %<>% subset(idents=opt$subset)
    }
    SeuratObject::Idents(srt) <- opt$celltype
    celltypes <- SeuratObject::Idents(srt) %>% levels()
    if(!is.null(opt$used_celltype)){
        used_celltype <- split(opt$used_celltype,',') %>% unlist()
        SeuratObject::Idents(srt) <- opt$celltype
        srt %<>% subset(idents=used_celltype)
        SeuratObject::Idents(srt) <- "CellType"
    }
    if(ncol(srt) > 2e+5){
        meta=srt@meta.data
        meta %<>% dplyr::slice_sample(prop=0.3,by=opt$celltype)
        srt %<>% subset(cells=rownames(meta))
    }
    fdata <- data.frame(gene_short_name=rownames(srt),row.names = rownames(srt))
    cds <- monocle::newCellDataSet(
        cellData=SeuratObject::GetAssayData(srt,'RNA'),
        phenoData = new("AnnotatedDataFrame", data = srt@meta.data),
        featureData = new("AnnotatedDataFrame", data = fdata),
        expressionFamily = VGAM::negbinomial.size()
    )
    hvg <- srt@assays$SCT@var.features
}

cds %<>% BiocGenerics::estimateSizeFactors() %>% BiocGenerics::estimateDispersions()

## 使用monocle选择的高变基因
# hvg <- monocle::dispersionTable(cds) %>%
#     dplyr::filter(mean_expression >= 0.1, dispersion_empirical >= 1 * dispersion_fit) %>%
#     dplyr::pull("gene_id")
# cds %<>% monocle::setOrderingFilter(hvg)
# pdf("plot_ordering_genes.pdf")
# plot_ordering_genes(cds)
# dev.off()

# deg_table <- monocle::differentialGeneTest(cds, fullModelFormulaStr = paste0("~",opt$celltype), cores = 1)
# print(head(deg_table))

## 差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
# deg_table %<>% 
#     dplyr::filter(qval < 0.01) %>%
#     dplyr::arrange(qval)

## 差异基因的结果文件保存
# deg_table %>%
#     dplyr::select(gene_short_name, everything()) %>%
#     data.table::fwrite(file = "monocle_DEG.xls",sep='\t')

## 轨迹构建基因可视化
cds %<>% monocle::setOrderingFilter(hvg)
# 这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
# setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
pdf("train.ordergenes.pdf")
monocle::plot_ordering_genes(cds)
dev.off()
# 出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)
# ordergene <- hvg[1:400] %>% na.omit()
# message(stringr::str_glue("Number of DEGs: {length(ordergene)}"))

# Step 2: 降维
cds %<>% monocle::reduceDimension(max_components = 2, method = "DDRTree")
# Step 3: 拟时间轴轨迹构建和在拟时间内排列细胞
cds %<>% monocle::orderCells(root_state = NULL)

# install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz',repos=NULL)
p_pseudotime <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE) +
    ggplot2::scale_colour_gradient2(low=hue[1],high=hue[4])
gs(p_pseudotime,name='Trajectory_Pseudotime',outdir=outdir,w=7,h=6)

# p <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", markers="NUPR1", size = 1)
# gs(p, name='Trajectory_Pseudotime-NUPR1',outdir=outdir,w=7,h=6)

p_cluster <- monocle::plot_cell_trajectory(cds, color_by = "CellType", size = 1, show_backbone = TRUE) +
    ggplot2::scale_colour_manual(values=hue)
gs(p_cluster,name='Trajectory_Cluster',outdir=outdir,w=7,h=6)

p_state <- monocle::plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
    ggplot2::scale_colour_manual(values=hue)
gs(p_state,name='Trajectory_State',outdir=outdir,w=7,h=6)

p_state_facet <- monocle::plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1) +
    ggplot2::scale_colour_manual(values=hue)
gs(p_state_facet,name='Trajectory_State_facet',w=14,outdir=outdir)

# 树形图
p1 <- monocle::plot_cell_trajectory(cds, x = 1, y = 2, color_by = "CellType") +
    theme(legend.position = "none", panel.border = element_blank()) + # 去掉第一个的legend
    ggplot2::scale_colour_manual(values=hue)
p2 <- plot_complex_cell_trajectory(cds,
    x = 1, y = 2,
    color_by = "CellType"
) +
    ggplot2::scale_colour_manual(values=hue) +
    theme(legend.title = element_blank())
gs(p1 | p2,name='Tree',w=14,outdir=outdir)

# 时间轴的细胞密度图
cds_pdata <- Biobase::pData(cds)
p_desity <- ggplot(cds_pdata, aes(Pseudotime, colour = CellType, fill = CellType)) +
    geom_density(bw = 0.5, size = 1, alpha = 0.7) +
    ggpubr::theme_classic2()+
    ggplot2::scale_fill_manual(values=hue)
gs(p_desity,name='Pseudotime_Density',w=14,outdir=outdir)

# 提取感兴趣的细胞（进行后续分析）
# 比如对State7的细胞感兴趣
# s.cells <- subset(pdata, State == "1") %>% rownames()
# saveseurat(s.cells, file = "Monocle2_State_0.seurat")
# # 保存结果
# write.csv(pData(cds), "pseudotime.csv")

# 指定基因的可视化
## 选择前4个top基因并将其对象取出

keygenes <- head(hvg, 3)
keygenes <- c('NUPR1',head(hvg, 3))
cds_subset <- cds[keygenes, ]
## 可视化：以state/Cluster/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")+ggplot2::scale_colour_manual(values=hue)
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "CellType")+ggplot2::scale_colour_manual(values=hue)
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")+ggplot2::scale_color_gradient(low=hue[1],high=hue[3])
gs(p1 | p2 | p3,name='Genes_pseudotimeplot',w=16, h=8,outdir=outdir)


# # 指定基因
# s.genes <- c("SELL", "CCR7", "IL7R", "CD84", "CCL5", "S100A4")
# p1 <- plot_genes_jitter(cds[s.genes, ], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes, ], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes, ], color_by = "State")
# plotc <- p1 | p2 | p3
# ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)
# 拟时序展示单个基因表达量
# colnames(pData(cds))
# pData(cds)$CCL5 <- log2(exprs(cds)["CCL5", ] + 1)
# p1 <- plot_cell_trajectory(cds, color_by = "CCL5") + scale_color_gsea()
# pData(cds)$S100A4 <- log2(exprs(cds)["S100A4", ] + 1)
# p2 <- plot_cell_trajectory(cds, color_by = "S100A4") + scale_color_gsea()
# gs(p1 + p2,name='_Single_Genes_pseudotimeplot',w=16, h=8)


# 8. 寻找拟时相关的基因（拟时差异基因）
# 这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
# 如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[na.omit(hvg), ],
    cores = 1,
    fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# Time_diff <- Time_diff[, c(5, 2, 3, 4, 1, 6)] # 把gene放前面，也可以不改
fwrite(Time_diff, "Time_diff_all.xls",sep='\t')
# Time_diff <- fread('Time_diff_all.csv')

Time_genes <- Time_diff %>%
    dplyr::slice_min(order_by=qval,n=49,with_ties = TRUE) %>% 
    pull('gene_short_name')
Time_genes <- c('NUPR1',Time_genes)
p <- monocle::plot_pseudotime_heatmap(cds[Time_genes, ],hmcols=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100), num_clusters = 4, show_rownames = T, return_heatmap = T)
gs(p,name='Time_heatmap',w=5, h=9,outdir=outdir)

# 单细胞轨迹的“分支”分析
p <- monocle::plot_cell_trajectory(cds, color_by = "State")+ggplot2::scale_colour_manual(values=hue)
gs(p,name='Time_Trajectory',w=8, h=7,outdir=outdir)

base::saveRDS(cds, file = paste0(opt$filename,".rds"))