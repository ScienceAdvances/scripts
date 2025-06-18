scone <- function(
    object=srt,
    condaenv = "/home/zktbk0/.local/share/mamba/envs/R",
    hue = using::hue("NPG"),
    batch_key = "Sample",
    outdir,
    nfeatures = 3000,
    resolution = 0.8,
    vars2regress = NULL,
    point_size = NULL,
    integrate_method = Seurat::HarmonyIntegration, # CCAIntegration, JointPCAIntegration, RPCAIntegration
    npcs = 50,
    ndims = 20,
    w=9,
    h=8,
    verbose = TRUE) {
    using::using(using,reticulate,sctransform,patchwork)
    reticulate::use_condaenv(condaenv = condaenv, required = TRUE)

    object %<>% Seurat::NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
        Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
    Seurat::VariableFeaturePlot(
        object,
        cols = c("black", "red"),
        pt.size = 1,
        log = NULL,
        selection.method = NULL,
        assay = NULL,
        raster = FALSE
    ) %>% 
    using::gs(outdir = outdir, name = "VariableFeaturePlot", w = 8, h = 7)
    # object[["RNA"]] <- split(object[["RNA"]], f = dplyr::pull(object@meta.data, batch_key))
    # object %<>% Seurat::SCTransform(vars.to.regress = vars2regress, variable.features.n = nfeatures, vst.flavor = "v2", verbose = verbose)
    # SeuratObject::DefaultAssay(object) <- "SCT"
    # features <- SeuratObject::VariableFeatures(object, method = "sct")
    object %<>% Seurat::ScaleData(features = Seurat::VariableFeatures(object), vars.to.regress = vars2regress) %>%
        Seurat::RunPCA(npcs = npcs, verbose = verbose) %>%
        Seurat::RunUMAP(dims = 1:ndims, reduction = "pca", verbose = verbose) %>%
        Seurat::RunTSNE(dims = 1:ndims, reduction = "pca", verbose = verbose)
    # 查看批次效应
    p <- scCustomize::DimPlot_scCustom(object, reduction = "umap", label.size = 6, group.by = batch_key, label = F, figure_plot = TRUE, colors_use = hue) &
        ggplot2::ggtitle("Batch_Effect Before") & 
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4)) & 
        ggplot2::theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
    using::gs(p,outdir = outdir, name = "UMAP_Batch_Effect_before", w = w, h = 12)
    p <- scCustomize::DimPlot_scCustom(object, reduction = "tsne", label.size = 6, group.by = batch_key, label = F, figure_plot = TRUE, colors_use = hue) &
        ggplot2::ggtitle("Batch_Effect Before") & theme(plot.title = element_text(hjust = 0.5)) & 
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4)) & 
        ggplot2::theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
    using::gs(p,outdir = outdir, name = "TSNE_Batch_Effect_before", w = w, h = 12)

    object[["RNA"]] <- split(object[["RNA"]], f = dplyr::pull(object@meta.data, batch_key))
    # 去除批次效应
    after <- Seurat::IntegrateLayers(
        object = object,
        method = integrate_method,
        orig.reduction = "pca",
        new.reduction = "integrated",
        normalization.method = "LogNormalize", #"SCT",
        verbose = verbose
        # k.weight=35 # 如果报错需要修改
    )

    after %<>% 
        Seurat::RunUMAP(dims = 1:ndims, reduction = "integrated", verbose = verbose) %>%
        Seurat::RunTSNE(dims = 1:ndims, reduction = "integrated", verbose = verbose) %>%
        Seurat::FindNeighbors(dims = 1:ndims, reduction = "integrated", verbose = verbose) %>%
        Seurat::FindClusters(resolution = resolution,random.seed=1,  algorithm=4, verbose = verbose)
    SeuratObject::Idents(after) <- "seurat_clusters"

    p <- scCustomize::DimPlot_scCustom(after, reduction = "umap", label.size = 6, group.by = batch_key, label = F, figure_plot = TRUE, colors_use = hue) &
        ggplot2::ggtitle("Batch_Effect After") & 
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4)) & 
        ggplot2::theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
    using::gs(p,outdir = outdir, name = "UMAP_Batch_Effect_after", w = w, h = h)
    p <- scCustomize::DimPlot_scCustom(after, reduction = "tsne", label.size = 6, group.by = batch_key, label = F, figure_plot = TRUE, colors_use = hue) &
        ggplot2::ggtitle("Batch_Effect After") & theme(plot.title = element_text(hjust = 0.5)) & 
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 4)) & 
        ggplot2::theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
    using::gs(p,outdir = outdir, name = "TSNE_Batch_Effect_after", w = w, h = h)

    # for (i in c(base::seq(0.01, 0.1, 0.02), base::seq(0.1, 1, 0.1))) {
    #     after %<>% Seurat::FindClusters(resolution = i, cluster.name = base::paste("Cluster", i, sep = "_"), verbose = verbose)
    # }
    # clusters <- base::grep(x = colnames(after@meta.data), pattern = "^Cluster", value = TRUE)

    # p5 <- scCustomize::DimPlot_scCustom(after, reduction = "umap", label.size = 6, repel = T, group.by = clusters, split_seurat = T, label = T, figure_plot = TRUE, colors_use = hue) &
    #     Seurat::NoLegend()
    # using::gs(p5, outdir = outdir, name = "UMAP_Resolutions", w = 32, h = base::ceiling(length(clusters) / 3) * 6)

    # p6 <- scCustomize::DimPlot_scCustom(after, reduction = "tsne", label.size = 6, repel = T, group.by = clusters, split_seurat = T, label = T, figure_plot = TRUE, colors_use = hue) &
    #     Seurat::NoLegend()
    # using::gs(p6, outdir = outdir, name = "TSNE_Resolutions", w = 32, h = base::ceiling(length(clusters) / 3) * 6)

    # p7 <- clustree::clustree(after@meta.data, prefix = "Cluster_")
    # using::gs(p7, outdir = outdir, name = "Resolution_Tree", w = 16, h = 10)

    p <- scCustomize::DimPlot_scCustom(after, reduction = "umap", label.size = 6, repel = T, group.by = "seurat_clusters", label = T, figure_plot = TRUE, colors_use = hue)&
        Seurat::NoLegend()
    using::gs(p,outdir = outdir, name = "UMAP_Cluster", w = 8, h = 7)
    p=scCustomize::DimPlot_scCustom(after, reduction = "tsne", label.size = 6, repel = T, group.by = "seurat_clusters", label = T, figure_plot = TRUE, colors_use = hue) &
        Seurat::NoLegend()
    using::gs(p,outdir = outdir, name = "TSNE_Cluster", w = 8, h = 7)
    # 标准化
    SeuratObject::DefaultAssay(after) <- "RNA"
    #     Seurat::RunPCA(features = features) %>%
    #     Seurat::JackStraw(num.replicate = 100) %>%
    #     Seurat::ScoreJackStraw(dims = 1:ndims)
    
    reshape2::dcast(as.data.frame(table(data.frame("Cluster" = SeuratObject::Idents(after), "Sample" = after$Sample))), Sample ~ Cluster) %>%
        data.table::fwrite(file.path(outdir, "cluster_summary.xls"), sep = "\t")

    # PCA plot
    pElbowPlot <- Seurat::ElbowPlot(after)$data %>%
        ggplot2::ggplot(fill = hue) +
        ggplot2::geom_point(aes(x = dims, y = stdev)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "PC", y = "Standard Deviation") +
        ggplot2::theme(text = element_text(size = 18))
    using::gs(pElbowPlot, outdir = outdir, name = "PCA_ElbowPlot", w = 7, h = 7)

    # pJackStrawPlot <- Seurat::JackStrawPlot(after) +
    #     ggplot2::theme_bw() +
    #     ggplot2::scale_color_manual(values = hue) +
    #     ggplot2::theme(text = element_text(size = 18))
    # using::gs(pJackStrawPlot, outdir = outdir, name = "PCA_JackStrawPlot", w = 7, h = 7)
    after %<>% JoinLayers()
    return(after)
}

