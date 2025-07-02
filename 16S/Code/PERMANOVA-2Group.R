RareAbundance=mp_extract_assays(mpse, .abundance='RareAbundance', byRow = TRUE)
dist_matrix <- vegdist(t(RareAbundance), method = "bray")
group <- metadata$Group
# 执行PERMANOVA（使用adonis2函数）
permanova_res <- adonis2(dist_matrix ~ group, permutations = 9999, method = "bray")

# 提取统计结果
R_squared <- round(permanova_res$R2[1], 3)  # R²值
p_value <- round(permanova_res$`Pr(>F)`[1], 3)  # p值

# 生成所有可能的组间对比标签（示例分组为UC和对照组）
distance_detail <- dist_matrix %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Sample1") %>% 
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Distance") %>% 
  dplyr::filter(Sample1 < Sample2) %>% 
  dplyr::mutate(
    Group1 = str_remove(Sample1,'\\d{1,2}_T'),  # 假设UC为实验组
    Group2 = str_remove(Sample2,'\\d{1,2}_T'),
    # 生成完整的对比标签（包含组内和所有组间组合）
    Comparison = case_when(
      Group1 == "HC" & Group2 == "HC" ~ "HC",
      Group1 == "UC" & Group2 == "UC" ~ "UC",
      # Group1 == "UCPC" & Group2 == "UCPC" ~ "UCPC",
      TRUE ~ "Between groups"  # 其他情况均为组间
    )
  ) %>% 
  # 强制标签顺序（确保UC和Control组内在前，组间在后）
  dplyr::mutate(Comparison = factor(Comparison, levels = c("HC", "UC", "Between groups")))

distance_detail$Comparison %>% table()
# 绘制箱线图（带统计量注释）
p=ggplot(distance_detail, aes(x = Comparison, y = Distance, fill = Comparison)) +
  geom_boxplot(width = 0.6, outlier.size = 1) +
  scale_fill_manual(values = hue) +
  labs(
    x = "", 
    y = "Bray-Curtis Distance",
    title = "PERMANOVA: Group Distance Distribution",
    caption = paste0("PERMANOVA: R²=", R_squared, ", p=", p_value)
  ) +
  # 调整坐标轴和文本样式
  scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, 0.05)) +  # 匹配示例图y轴范围
  theme_bw(base_size = 18) + setheme()
using::gs(p, outdir = outdir, name = 'PERMANOVA', w = 7, h = 6)
