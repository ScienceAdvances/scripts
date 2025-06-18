marker_dotplot <- function(object, marker_anno, outdir, value_colors= 'RdBu', hue,w,h) {
    # https://mp.weixin.qq.com/s/2VJrpk5B0Q3yb65mnoMPzQ
    data <- marker_anno %>%
        dplyr::as_tibble() %>%
        tidyr::separate_rows(Marker, sep = ",") %>% 
        dplyr::distinct(Marker,.keep_all = TRUE) %>% 
        tidyr::drop_na(Marker)
    p1 <- Seurat::DotPlot(
        object,
        features = split(data$Marker, data$Anno),
        cols = value_colors
    ) +
        Seurat::RotatedAxis() + # 来自Seurat
        # scale_color_gradientn(colours = value_colors)+
        theme(
            panel.border = element_rect(color = "black"),
            panel.spacing = unit(1, "mm"),
            axis.title = element_blank(),
            axis.text.y = element_blank()
        )
    # hue <- c("#4A9D47", "#F19294", "#E45D61", "#96C3D8", "#5F9BBE", "#F5B375", "#C0937E", "#67A59B", "#A5D38F")
    df <- data.frame(x = 0, y = levels(object), stringsAsFactors = F)
    df$y <- factor(df$y, levels = df$y)
    p0 <- ggplot(df, aes(x, y, color = factor(y))) +
        geom_point(size = 6, show.legend = F) +
        scale_color_manual(values = hue) +
        theme_classic() +
        scale_x_continuous(expand = c(0, 0)) +
        theme(
            plot.margin = margin(r = 0),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.ticks = element_blank(),
            axis.line = element_blank()
        )

    p <- plot_grid(p0, p1, align = "h", axis = "bt", rel_widths = c(1.5, 8.5))+theme(plot.margin=margin(0.1, 0.1,0.1,0.6, unit = "cm"))#top, right, bottom, left
    gs(p,outdir=outdir,name='marker_dotplot',w=w,h=h)
}
