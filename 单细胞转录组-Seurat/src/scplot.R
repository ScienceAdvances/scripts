scplot <- function(object, bunch, hue, outdir, marker, assay = "RNA") {
  Seurat::DefaultAssay(object) <- assay
  nsample <- length(unique(object$Sample))
  ngroup <- length(unique(object$Group))
  ncluster <- object@meta.data %>%
    dplyr::pull(bunch) %>%
    unique() %>%
    length()
  p1 <- plot.cluster.std(data = object, clusters = bunch, xlab = "Cluster number", log = FALSE, group = "Group", legend.title = "Group", widths = c(2, 2), hue = hue)
  using::gs(p1,outdir=outdir,name="cluster_number",w = 10,h = 6)
  # 按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
  temp <- data.frame(table(object@meta.data %>% dplyr::pull(bunch), object$Sample))
  # 按样本/组展示细胞占比。柱状图。横坐标是样本/组，纵坐标是细胞占比
  colnames(temp) <- c("Cluster", "Sample", "number")
  p <- ggplot(temp, aes(x = Sample, fill = Cluster, y = number)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    labs(x = "", y = "") +
    guides(fill = guide_legend(title = "", reverse = TRUE)) +
    scale_fill_manual(values = hue) +
    theme_bw() +
    SCTheme+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  using::gs(p,outdir=outdir,name="cluster_sample_cellcounts",w = 12,h = 8)

  temp <- data.frame(table(object@meta.data %>% dplyr::pull(bunch), object$Group))
  colnames(temp) <- c("Cluster", "Group", "number")
  p <- ggplot(temp, aes(x = Group, fill = Cluster, y = number)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    labs(x = "", y = "") +
    guides(fill = guide_legend(title = "", reverse = TRUE)) +
    scale_fill_manual(values = hue) +
    theme_bw() + SCTheme
  using::gs(p,outdir=outdir,name="cluster_group_cellcounts",w = 12,h = 8)
  # 分样本展示
  p <- scCustomize::DimPlot_scCustom(object,reduction = 'umap',group.by = 'Sample', label = FALSE,colors_use=hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_umap.sample',w=12,h=8)
  p <- scCustomize::DimPlot_scCustom(object,reduction = 'tsne',group.by = 'Sample', label = FALSE,colors_use=hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_tsne.sample',w=12,h=8)

  # 分组展示
  scCustomize::DimPlot_scCustom(object,reduction = 'umap',group.by = 'Group', label = FALSE,colors_use=hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_tsne.group',w=12,h=8)
  scCustomize::DimPlot_scCustom(object,reduction = 'tsne',group.by = 'Group', label = FALSE,colors_use=hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_tsne.group',w=12,h=8)

  # CellType
  p <- scCustomize::DimPlot_scCustom(object, reduction = "umap",group.by = bunch, label = FALSE, colors_use = hue,figure_plot=TRUE) 
  using::gs(p,outdir=outdir,name='cluster_umap',w=12,h=8)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "tsne",group.by = bunch, label = FALSE, colors_use = hue) 
  using::gs(p,outdir=outdir,name='cluster_tsne',w=12,h=8)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "umap",split_seurat = T,group.by = bunch, label = T, colors_use = hue,figure_plot=TRUE) & Seurat::NoLegend() 
  using::gs(p,outdir=outdir,name='cluster_umap_label',w=12,h=8)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "tsne",group.by = bunch, label = T, colors_use = hue) & Seurat::NoLegend()
  using::gs(p,outdir=outdir,name='cluster_tsne_label',w=12,h=8)

  # 分样本展示 按样本分页;分组展示，按组分页
  if(ngroup>5){
    sw <- 25
    scolumns <- 5
    sh <- ceiling(nsample / 5) * 8
  }else{
    sw <- nsample*5
    scolumns <- nsample
    sh <- 6
  }
  if(ngroup>5){
    gw <- 25
    gcolumns <- 5
    gh <- ceiling(ngroup / 5) * 8
  }else{
    gw <- ngroup*5
    gcolumns <- ngroup
    gh <- 6
  }
  p <- scCustomize::DimPlot_scCustom(object, reduction = "umap", split.by = "Sample", label = FALSE, split_seurat=TRUE,num_columns = scolumns, colors_use = hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_umap_splitsample',w=sw,h=sh)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "tsne", split.by = "Sample", label = FALSE, split_seurat=TRUE,num_columns = scolumns, colors_use = hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_tsne_splitsample',w=sw,h=sh)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "umap", split.by = "Group", label = FALSE, split_seurat=TRUE,num_columns = gcolumns, colors_use = hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_umap_splitgroup',w=gw,h=gh)
  p <- scCustomize::DimPlot_scCustom(object, reduction = "tsne", split.by = "Group", label = FALSE, split_seurat=TRUE,num_columns = gcolumns, colors_use = hue,figure_plot=TRUE)
  using::gs(p,outdir=outdir,name='cluster_tsne_splitgroup',w=gw,h=gh)

  # marke基因展示
  marker_number <- table(marker$Cluster) %>% reshape2::melt()
  marker_number$Cluster <- as.factor(marker_number$Var1)
  p <- ggplot(data = marker_number, mapping = aes(x = Cluster, y = value, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.8) +
    theme_classic() +
    ylab("Number") +
    scale_fill_manual(values = hue) +
    geom_text(aes(label = value), size = 3, hjust = 0.5, vjust = -0.5, position = "stack")
  using::gs(p, outdir=outdir,name='marker_number',w=14,h=6)

  # top基因展示
  all_top10_markers <- marker %>%
    dplyr::filter(!is.infinite(avg_log2FC)) %>%
    dplyr::group_by(Cluster) %>%
    top_n(n = 1, wt = avg_log2FC) %>%
    as.data.frame() %>%
    dplyr::distinct(., Feature, .keep_all = T) %>%
    head(n = 10)
  height <- ceiling(length(all_top10_markers$Feature) / 5) * 4
  p <- scCustomize::VlnPlot_scCustom(object, features = all_top10_markers$Feature, num_columns=5, colors_use=hue, plot_boxplot = TRUE)
  using::gs(p,outdir=outdir,name="top10gene_vilion",w = 25,h = 7)

  p <- scCustomize::DotPlot_scCustom(object, features = unique(unique(all_top10_markers$Feature)),colors_use=dot_hue, flip_axes = T, x_lab_rotate = TRUE)
    using::gs(p,outdir=outdir,name="cluster_top1_dotplot",w = 12,h = 12)

  p <- scCustomize::FeaturePlot_scCustom(object, features = all_top10_markers$Feature, min.cutoff = "q9",num_columns=5,figure_plot=T,order=T,colors_use=c("lightgrey", "red"))
    using::gs(p,outdir=outdir,name="cluster_top1_umap",w = 18,h = height)

  cluster_top10_markers <- marker %>%
    dplyr::group_by(Cluster) %>%
    dplyr::top_n(n = 6, wt = avg_log2FC)
  features <- unique(cluster_top10_markers$Feature)
  tmp <- subset(object, downsample = 2000)
  SeuratObject::DefaultAssay(tmp) <- assay
  tmp %<>% Seurat::ScaleData(features = features, assay=assay)
  p <- Seurat::DoHeatmap(tmp, features = features, group.colors = hue) + NoLegend() + theme(axis.text.y = element_text(size = 6))+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))
  using::gs(p, outdir=outdir,name="cluster_top6_markers_heatmap",w = ceiling(ncluster / 20) * 8,h = 10)
}



plot.clusters.group <- function(data, clusters,
                                group, widths = c(3, 1), log = TRUE,
                                legend.title = "Group", hue, xlab = "") {
    ## take an integrated Seurat object, plot distributions over orig.ident
    mytheme <- theme(
        plot.title = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
    )

    count_table <- table(data@meta.data[, clusters], data@meta.data[, group])
    count_mtx <- as.data.frame.matrix(count_table)
    count_mtx$Cluster <- rownames(count_mtx)
    melt_mtx <- melt(count_mtx)
    melt_mtx$Cluster <- as.factor(melt_mtx$Cluster)
    cluster_size <- aggregate(value ~ Cluster, data = melt_mtx, FUN = sum)

    if ("0" %in% cluster_size$Cluster) {
        sorted_labels <- paste(cluster_size$Cluster[order(cluster_size$value)])
    } else {
        sorted_labels <- paste(cluster_size$Cluster[order(cluster_size$value)])
    }

    cluster_size$Cluster <- factor(cluster_size$Cluster, levels = sorted_labels)
    melt_mtx$Cluster <- factor(melt_mtx$Cluster, levels = sorted_labels)

    colnames(melt_mtx)[2] <- "dataset"

    ################### p1
    head(cluster_size)
    if (log) {
        p1 <- ggplot(cluster_size, aes(y = Cluster, x = value)) +
            geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
            theme_bw() +
            scale_x_log10() +
            xlab("Cells per Cluster, log10 scale") +
            ylab("") +
            mytheme
    } else {
        p1 <- ggplot(cluster_size, aes(y = cluster, x = value)) +
            geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
            theme_bw() +
            xlab("Cells per Cluster") +
            ylab("") +
            mytheme
    }
    ################### p2
    ########################### color 1
    p2 <- ggplot(melt_mtx, aes(x = Cluster, y = value, fill = dataset)) +
        geom_bar(position = "fill", stat = "identity", ) +
        theme_bw() +
        coord_flip() +
        scale_fill_manual(values = hue) +
        ylab(paste0("Fraction of cells in each ", tolower(legend.title))) +
        xlab(xlab) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent, expand = c(0.01, 0.01)) +
        mytheme
    p3 <- ggplot(melt_mtx, aes(x = Cluster, y = value, fill = dataset)) +
        geom_bar(stat = "identity", ) +
        theme_bw() +
        coord_flip() +
        scale_fill_manual(values = hue) +
        ylab(paste0("count of cells in each ", tolower(legend.title))) +
        xlab(xlab) +
        theme(legend.position = "top") +
        guides(fill = guide_legend(title = legend.title)) +
        mytheme
    wrap_plots(ncol = 2, p2, p1, p3, widths = widths)
}

plot.cluster.std <- function(data, clusters, group = orig.ident, widths = c(2, 2), log = TRUE, legend.title = "Group", hue, xlab = "") 
{
    mytheme <- theme(
        plot.title = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        # axis.line = element_line(color = "black"),
        # axis.ticks = element_line(color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid=element_blank(), # 去网格线
        # legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        # axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
    )

    count_table <- table(data@meta.data[, clusters], data@meta.data[, group])
    count_mtx <- as.data.frame.matrix(count_table)
    count_mtx$Cluster <- rownames(count_mtx)
    melt_mtx <- melt(count_mtx)
    melt_mtx$Cluster <- as.factor(melt_mtx$Cluster)
    cluster_size <- aggregate(value ~ Cluster, data = melt_mtx, FUN = sum)

    if ("0" %in% cluster_size$Cluster) {
        # sorted_labels <- paste(sort(as.integer(levels(cluster_size$Cluster)),decreasing = T))
        sorted_labels <- paste(cluster_size$Cluster[order(cluster_size$value)])
    } else {
        sorted_labels <- paste(cluster_size$Cluster[order(cluster_size$value)])
    }

    cluster_size$Cluster <- factor(cluster_size$Cluster, levels = sorted_labels)
    melt_mtx$Cluster <- factor(melt_mtx$Cluster, levels = sorted_labels)

    colnames(melt_mtx)[2] <- "dataset"

    ################### p1
    head(cluster_size)
    if (log) {
        p1 <- ggplot(cluster_size, aes(y = Cluster, x = value)) +
            geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
            theme_bw() +
            scale_x_log10() +
            xlab("Cells per cluster, log10 scale") +
            ylab("") +
            mytheme +
            scale_fill_manual(values = hue)
    } else {
        p1 <- ggplot(cluster_size, aes(y = Cluster, x = value)) +
            geom_bar(position = "dodge", stat = "identity", fill = "grey60") +
            theme_bw() +
            xlab("Cells per Cluster") +
            ylab("") +
            mytheme +
            scale_fill_manual(values = hue)
    }

    p2 <- ggplot(melt_mtx, aes(x = Cluster, y = value, fill = dataset)) +
        geom_bar(position = "fill", stat = "identity", ) +
        theme_bw() +
        coord_flip() +
        scale_fill_brewer(palette = "Set2") +
        ylab("Proportion of cells") +
        xlab(xlab) +
        theme(legend.position = "left") +
        guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent, expand = c(0.01, 0.01)) +
        mytheme +
        scale_fill_manual(values = hue)
    wrap_plots(ncol = 2, p2, p1, widths = widths)
}
SCTheme <- theme(
    axis.title = element_text(
      face = "bold",
      size = "16", color = "black"
    ),
    # legend.position = 'right',
    axis.text.x = element_text(
      color = "black", face = "bold",
      size = 10, hjust = 0.5, vjust = 0.5
    ),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(face = "bold", color = "black", size = 10),
    legend.title = element_text(face = "bold", color = "black", size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "right",
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(colour = "white", fill = "grey"),
    plot.title = element_text(face = "bold", color = "black", lineheight = .8, hjust = 0.5, size = 5)
)
