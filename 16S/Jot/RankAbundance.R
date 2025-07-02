otu_relative <- t(apply(otu, 1, function(x) x / sum(x)))
for (i in colnames(otu)) {
  otu[, i] <- rev(sort(otu_relative[, i]))
}
library(ggplot2)
library(reshape2)

# 转换数据格式
otu_long <- reshape2::melt(otu, id.vars = "OTU")
colnames(otu_long) <- c('Sample','Feature','value')
# 添加排名列
otu_long$Rank <- ave(otu_long$value, otu_long$Sample, 
                     FUN = function(x) rank(-x, ties.method = "first"))
df <- merge(otu_long, metadata, by.x = 'Sample', by.y = 0)

otu_long %>% head(2)
# 绘制曲线（分样本颜色）
p <- ggplot(df, aes(x = Rank, y = value, color = Group)) +
  geom_line(aes(group=Sample),linewidth = 0.8, alpha = 0.7) +
  scale_y_log10(labels = scales::percent) +  # 纵轴对数转换
  scale_color_manual(values = hue4[2:1])+
  labs(x = "Species Rank", y = "Relative Abundance (%)") +
  theme_bw() +
  setheme()
using::gs(p, outdir = outdir, name = "RankAbundance", w = 12, h = 6)
