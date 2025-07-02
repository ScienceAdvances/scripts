library(tidyverse)
source("Code/src/RConfig.R")
outdir <- "Result/02.稀释曲线"
using::mkdir(outdir)

# 运行计算
shannon_curves <- alpha_curves(otu, step = 1000, method = "shannon")
df <- merge(shannon_curves, metadata, by.x = 'Sample', by.y = 0)

# 计算均值±标准差
stats <- df %>%
  dplyr::group_by(Group, Depth) %>%
  dplyr::summarise(Mean = mean(Alpha), SD = sd(Alpha))

# 绘制带置信区间的曲线
p <- ggplot(stats, aes(x = Depth, y = Mean, color = Group)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = hue4[2:1]) +
  scale_fill_manual(values = hue4[2:1]) +
scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())  # 自动添加K/M/G单位[3](@ref)
  )+
  labs(x = "Sequencing Depth", y = "Shannon Diversity Index") +
  theme_bw(base_size = 12) + setheme(45)
using::gs(p, outdir = outdir, name = "ShannonDiversityIndex-Group", w = 7, h = 6)


# 计算均值±标准差
stats <- df %>%
  group_by(Sample, Depth) %>%
  dplyr::summarise(Mean = mean(Alpha),
            SD = sd(Alpha)) %>% dplyr::mutate(Depth= as.numeric(Depth))

# 绘制带置信区间的曲线
p <- ggplot(stats, aes(x = Depth, y = Mean, color = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Sample), 
              alpha = 0.2, color = NA) +
  scale_fill_manual(values = hue4) +
  scale_color_manual(values = hue4) +
scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())  # 自动添加K/M/G单位[3](@ref)
  )+
  labs(x = "Sequencing Depth", y = "Shannon Diversity Index") +
  theme_bw(base_size = 12) + setheme(45)
using::gs(p, outdir = outdir, name = "ShannonDiversityIndex-Sample", w = 7, h = 6)
