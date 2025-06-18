source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(data.table,tidyverse)
hue_day <- c('#F6E0D2','#DFA398','#9C6755','#659794','#EA967C','#F5C98E','#D65B5A','#586085')

fdata <- data.table::fread("TrinityRaw/corset/counts.txt") %>% tibble::column_to_rownames("V1")
pdata <- data.table::fread("Meta/pdata.tsv") %>% tibble::column_to_rownames("Sample")
# A1命名：CK5-6
# A2命名：CK7-8
# A3命名：CK9-10
# B1命名：NPA5-6
# B2命名：NPA7-8
# B3命名：NPA9-10

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01)


# 创建数据框
data <- data.frame(
  group = c('CK5-6vsNPA5-6', "CK7-8vsNPA7-8", "CK9-10vsNPA9-10"),
  up_regulated = c(nrow(g1[log2FoldChange > 2,]), nrow(g2[log2FoldChange > 2,]), nrow(g3[log2FoldChange > 2,])),
  down_regulated = c(-nrow(g1[log2FoldChange < -2,]), -nrow(g2[log2FoldChange < -2,]), -nrow(g3[log2FoldChange < -2,]))
)

# 将数据转换为长格式，以便ggplot2使用
library(reshape2)
data_long <- reshape2::melt(data, id.vars = "group")

# 使用ggplot2绘制柱状图
p <- ggplot(data_long, aes(x = group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = hue_day,
                    labels = c("up - regulated", "down - regulated"),
                    name = "") +
  xlab("") +
  ylab("Numbers of DEGs") +
  theme_minimal()+Theme
using::gs(p,outdir='F2',name='A')

Theme <- theme(
    axis.title = element_text(
      face = "bold",
      size = "18", color = "black"
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

