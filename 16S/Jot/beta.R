source("Code/RConfig.R")
outdir <- "Result/Beta"
using::mkdir(outdir)

data$cal_betadiv()
t1 <- trans_beta$new(dataset = data, group = "Group", measure = "bray")
t1$cal_ordination(method = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
p1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))+
theme_bw()+theme4
using::gs(p1, outdir = outdir, name = 'PCoA', w = 8, h = 7)

t1$cal_ordination(method = "PCA")

t1 <- trans_beta$new(dataset = data, group = "Group", measure = "bray")
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
conflicts_prefer(dplyr::mutate)
p <- t1$plot_group_distance(add = "mean")+theme_bw()+theme4
using::gs(p, outdir = outdir, name = 'bray_boxplot', w = 8, h = 7)


# 3组
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "wilcox")
# parameters in plot_group_distance function will be passed to the plot_alpha function of trans_alpha class
t1$plot_group_distance(plot_type = "ggviolin", add = "mean_se")
p <- t1$plot_group_distance(add = "mean")+theme_bw()+theme4
using::gs(p, outdir = outdir, name = 'bray_boxplot_vs', w = 8, h = 7)
