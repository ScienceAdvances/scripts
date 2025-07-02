relative_abundance <- function(table){
    tibble::as_tibble(table,rownames='ASV') %>% 
        tidyr::pivot_longer(cols=-ASV) %>% 
        dplyr::group_by(name) %>% 
        dplyr::mutate(total=sum(value)) %>% 
        dplyr::mutate(relative_abun=value/total) %>% 
        tidyr::pivot_wider(id_cols=c("ASV"),names_from = name,values_from = relative_abun) %>% 
        dplyr::mutate_if(is.numeric,round)
}

betaplot <- function(obejct,ordination,distance,hue,w,h) {
    p <- obejct %>% mp_plot_dist(.distmethod = !!distance, .group = Group, group.test=TRUE, textsize=4) +
        scale_fill_manual(values=hue) +
        scale_color_manual(values=hue) + theme_bw() + setheme()
    using::gs(p, outdir = outdir, name = glue::glue('BetaDiversity-{distance}-boxplot'), w = 7, h = 6)
    p1 <- obejct %>%
            mp_plot_ord(
            .ord = !!ordination, 
            .group = Group, 
            .color = Group, 
            .size = Observe,
            .alpha = Observe,
            ellipse = TRUE,
            show.legend = FALSE # don't display the legend of stat_ellipse
            ) +
            scale_fill_manual(values=hue) +
            scale_color_manual(values=hue) 
    p2 <- obejct %>% 
            mp_plot_ord(
            .ord = !!ordination, 
            .group = Group, 
            .color = Group, 
            .size = Shannon, 
            .alpha = Shannon,
            ellipse = TRUE,
            show.legend = FALSE # don't display the legend of stat_ellipse 
            ) +
            scale_fill_manual(
            values = hue
            ) +
            scale_color_manual(
            values=hue
            )
    p <- p1 + p2 & theme_bw() & setheme()
    using::gs(p, outdir = outdir, name = glue::glue('BetaDiversity-{distance}-{ordination}'), w = w, h = h)
}

