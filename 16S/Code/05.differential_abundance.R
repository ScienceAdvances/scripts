source('Code/src/RConfig.R')
outdir <- "Result/05.differential_abundance"
using::mkdir(outdir)
META=Sys.getenv("META")
t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 8)
p <- t1$plot_bar(others_color = "grey70",color_values = hue200,facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE) + 
    theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "bar_Phylum", w = 12, h = 8)

t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 8)
p <- t1$plot_bar(bar_type = "notfull", 
    use_alluvium = TRUE, 
    clustering = TRUE, 
    xtext_angle = 30, 
    xtext_size = 3, 
    color_values = hue200)+
    theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "bar_Genus", w = 12, h = 8)

# stack plot
t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
p <- t1$plot_bar(others_color = "grey70",color_values = hue200, legend_text_italic = FALSE)+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "stackbar_Phylum", w = 6, h = 8)

# boxplot
t1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 15)
p <- t1$plot_box(group = "Group", xtext_angle = 30,color_values = hue200)+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "boxplot", w = 9, h = 6)
#  c("Phylum", "Class", "Order", "Family", "Genus")
# heatmap
t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 40)
p <- t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "heatmap", w = 12, h = 9)

t1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 40)
p <- t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "heatmap-Class", w = 12, h = 9)
t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 40)
p <- t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "heatmap-Genus", w = 12, h = 9)

# plot_line
t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 5, groupmean = "Group")
p <- t1$plot_line(color_values = hue200) + theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "line_Phylum", w = 9, h = 6)

t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 5, groupmean = "Group")
p <- t1$plot_line(position = position_dodge(0.3), xtext_angle = 0,color_values = hue200) + theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "line_Genus", w = 9, h = 6)

# pie
t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
p <- t1$plot_pie(facet_nrow = 1, add_label = TRUE,color_values = hue200)
using::gs(p, outdir = outdir, name = "pie_Phylum", w = 9, h = 6)

t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 8, groupmean = "Group")
p <- t1$plot_donut(label = TRUE,color_values = hue200) + theme_bw()
using::gs(p, outdir = outdir, name = "pie_Genus", w = 9, h = 6)

# radar
t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
p <- t1$plot_radar(values.radar = c("0%", "25%", "50%"),color_values = hue200, grid.min = 0, grid.mid = 0.25, grid.max = 0.5) + theme_bw()
using::gs(p, outdir = outdir, name = "radar_Phylum", w = 9, h = 6)

t1 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 10, groupmean = "Group")
p <- t1$plot_radar(values.radar = c("0%", "25%", "50%"),color_values = hue200, grid.min = 0, grid.mid = 0.25, grid.max = 0.5) + theme_bw()
using::gs(p, outdir = outdir, name = "radar_Genus", w = 9, h = 6)

# venn
dataset1 <- dataset$merge_samples("Group")
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
p <- t1$plot_venn() + ggplot2::scale_colour_manual(values=hue200)
using::gs(p, outdir = outdir, name = "venn", w = 7, h = 6)

# install.packages('ggtern')
# petal
if (unique(metadata$Group) %>% length() > 3) {
    dataset1 <- dataset$merge_samples("Group")
    t1 <- trans_venn$new(dataset1)
    t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8,petal_color=hue200, petal_color_center = "#BEBADA") + theme_bw()
    using::gs(p, outdir = outdir, name = "petal", w = 7, h = 6)
}

# tern
if (unique(metadata$Group) %>% length() == 3) {
    t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
    p <- t1$plot_tern()
    using::gs(p, outdir = outdir, name = "tern", w = 7, h = 6)
}


# require ggnested package; see https://chiliubio.github.io/microeco_tutorial/intro.html#dependence
test1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 10, high_level = "Phylum", delete_taxonomy_prefix = FALSE)
p1 <- test1$plot_bar(ggnested = TRUE, facet = 'Group', xtext_angle = 45)
using::gs(p1, outdir = outdir, name = "ggnested_bar", w = 16, h = 9)

# fixed subclass number in each phylum
test1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 30, show = 0, high_level = "Phylum", high_level_fix_nsub = 4)
p2 <- test1$plot_bar(ggnested = TRUE, xtext_angle = 30, facet = 'Group') + theme_bw() + setheme(xtext_rotation=45)
using::gs(p2, outdir = outdir, name = "fixed_bar", w = 16, h = 9)

# sum others in each phylum
test1 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 3, delete_taxonomy_prefix = FALSE)
p4 <- test1$plot_bar(ggnested = TRUE, high_level_add_other = TRUE, xtext_angle = 30, facet = 'Group')
using::gs(p4, outdir = outdir, name = "sum_bar", w = 16, h = 9)

t1 <- trans_abund$new(dataset = mt, ntaxa = 10, groupmean = "Group")
p <- t1$plot_bar(others_color = "grey70",color_values = hue200, legend_text_italic = FALSE)+
theme_bw() + setheme(xtext_rotation=45)
using::gs(p, outdir = outdir, name = "stackbar-Group", w = 6, h = 8)

dfs <- map(c("Phylum", "Class", "Order", "Family", "Genus"),function(x){
    d=dataset$taxa_abund[[x]]
    d2 <- apply(d > 0, 2, sum)
    d3 <- data.frame(Sample=names(d2),AAA = d2/sum(d2))
    colnames(d3) <- c('Sample',x)
    return(d3)
})
df <- purrr::reduce(dfs,inner_join,by='Sample')

df_long <- df %>% 
  pivot_longer(
    cols = -Sample,  # 保留Sample列，其他列转为长格式
    names_to = "Taxonomic_Level",  # 新列：分类层级（Phylum, Class...）
    values_to = "Count"  # 新列：对应数量
  )

# 将分类层级设为因子，固定顺序（避免按字母排序）
df_long$Taxonomic_Level <- factor(
  df_long$Taxonomic_Level,
  levels = c("Phylum", "Class", "Order", "Family", "Genus")
)
p <- ggplot(df_long, aes(x = Sample, y = Count, fill = Taxonomic_Level)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  ggplot2::scale_fill_manual(values=hue200) +  
  # 调整主题（参考网页8的简洁风格）
  theme_bw() + setheme(45)+
  labs(
    x = "",
    y = "Sequence Number Percent"
  ) 
using::gs(p, outdir = outdir, name = "AllClass", w = 16, h = 9)

p1 <- mpse %>% MicrobiotaProcess::mp_plot_dist(.distmethod = bray, .group = Group) %>% 
set_scale_theme(
      x = scale_fill_manual(
      values=hue, 
      guide = guide_legend(
      keywidth = 1, 
      keyheight = 0.5, 
      title.theme = element_text(size=8),
      label.theme = element_text(size=6)
   )
   ), 
      aes_var = Group # specific the name of variable 
   ) %>%
set_scale_theme(
   x = scale_color_gradient(
   guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
   ),
   aes_var = bray
   ) %>%
   set_scale_theme(
   x = scale_size_continuous(
   range = c(0.5, 3),
   guide = guide_legend(keywidth = 0.5, keyheight = 0.5)
   ),
   aes_var = bray
)
using::gs(p1, outdir = outdir, name = 'Heatmap', w = 9, h = 8)

# UPGMA
mpse %<>%
       mp_cal_clust(
         .abundance = 'RareAbundance', 
         distmethod = "bray",
         hclustmethod = "average", # (UPGMA)
         action = "add" # action is used to control which result will be returned
       )
sample.clust <- mpse %>% mp_extract_internal_attr(name='SampleClust')
p <- ggtree::ggtree(sample.clust) + 
       ggtree::geom_tippoint(aes(color=Group)) +
       scale_color_manual(values=hue) +
       ggtree::geom_tiplab(as_ylab = TRUE) +
       ggplot2::scale_x_continuous(expand=c(0, 0.01))
using::gs(p, outdir = outdir, name = 'UMPGMA', w = 7, h = 8)


