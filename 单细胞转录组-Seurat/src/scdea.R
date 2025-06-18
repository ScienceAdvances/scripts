scdea <- function(object, outdir,hue) {
    markers <- data.frame(Feature = "", p_val = "", avg_log2FC = "", pct.1 = "", pct.2 = "", p_val_adj = "", Cluster = "")[-1, ]
    for (i in levels(Idents(object))) {
        i_cluster_dir <- file.path(outdir, paste("cluster", i, sep = "_"))
        using::mkdir(i_cluster_dir)
        future::plan("multisession", workers = 4)
        cluster_markers <- Seurat::FindMarkers(object, ident.1 = i, verbose = FALSE)
        cluster_markers %<>% data.table::as.data.table(keep.rownames = "Feature")
        cluster_markers$Cluster <- i
        markers <- base::rbind(markers, cluster_markers)
        # cluster_dir_enrich=paste(cluster_dir,"enrichment",sep="/")
        # if(!file.exists(cluster_dir_enrich)){dir.create(cluster_dir_enrich)}
        if (nrow(cluster_markers) > 1) {
            # genelist=cluster_markers$gene
            # tmp = cluster_markers %>% dplyr::filter(p_val<0.05, avg_log2FC>0)
            # genelist=tmp$gene
            # try(enrichment(species=opt$type,outDir=cluster_dir_enrich,geneList=genelist))
            cluster_markers %>%
                dplyr::arrange(dplyr::desc(avg_log2FC), pct.1, p_val_adj) %>%
                data.table::fwrite(file.path(i_cluster_dir, "markers.xls"), sep = "\t")
            top20_markers <- cluster_markers %>%
                dplyr::top_n(n = 20, wt = avg_log2FC) %>%
                dplyr::pull("Feature")
            p <- scCustomize::VlnPlot_scCustom(object, features = top20_markers, num_columns = 5, colors_use = hue, plot_boxplot = TRUE)
            using::gs(p, outdir = i_cluster_dir, name = "top_gene_vilion", w = 25, h = 14)

            scCustomize::FeaturePlot_scCustom(object, features = top20_markers, min.cutoff = "q9", num_columns = 5, figure_plot = T, order = T, colors_use = c("lightgrey", "red")) %>%
                using::gs(outdir = i_cluster_dir, name = "top_gene_umap", w = 20, h = 14)

            scCustomize::DotPlot_scCustom(object, features = rev(unique(top20_markers)), flip_axes = T, x_lab_rotate = TRUE,colors_use = dot_hue) %>%
                using::gs(outdir = i_cluster_dir, name = "top_gene_dotplot", w = 15, h = 15)
        }
    }
    markers %>%
        dplyr::arrange(dplyr::desc(avg_log2FC), pct.1, p_val_adj) %>%
        data.table::fwrite(file.path(outdir, "all_markers.xls.gz"), sep = "\t")
    return(markers)
}
