source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(data.table,tidyverse)
hue_day <- c('#F6E0D2','#DFA398','#9C6755','#659794','#EA967C','#F5C98E','#D65B5A','#586085')

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)

genes <- c(g1,g2,g3) %>% unique()

d=fread('F2/数据/tf_classification.txt',header=FALSE)
# d$Gene <- str_remove(d$V1,'-[0-9]{1}[FR]{1}')
d$Gene <- map_chr(str_split(d$V1,'_'),~paste0(.x[1:4],collapse = '_'))

d1 <- d %>% dplyr::filter(Gene %in% genes)
d3 <- d1 %>% dplyr::group_by(V2) %>% dplyr::summarise(N=n()) %>% dplyr::arrange(desc(N)) %>% dplyr::filter(N>10)
d3$V2 <- factor(d3$V2,levels = d3$V2,ordered = T)
NPG=c("#E54B34","#4CBAD4","#009F86","#3B5387","#F29A7F","#8491B3","#91D1C1","#DC0000","#7E6047","#CCCCCC","#BC8B83","#33ADAD","#347988","#9F7685","#C1969A","#8BB0BB","#CE8662","#B04929","#A59487","#E3907E","#D46F5B","#41B4C1","#278C87","#726486","#DA988C","#88A0B7","#B9AC91","#C63517","#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#8EC0BE","#D85935","#985738","#B9AFA9","#E67059","#E5BFB9","#B2CED4","#779F99","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#7E7770","#CCCCCC")

p <- ggplot(d3, aes(x = V2, y = N, fill = V2)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(x = "Treatment", y = "Number of TFs and TRs", fill = "TF/TR type") +
  theme_minimal()+Theme+xlab('')+ggplot2::scale_fill_manual(values = c(hue_day,NPG))
using::gs(p,outdir='F2',name='D',w=16,h=9)

Theme <- theme(
    axis.title = element_text(
      face = "bold",
      size = "18", color = "black"
    ),
    # legend.position = 'right',
    axis.text.x = element_text(
      color = "black", face = "bold",
      size = 10, 
    #   hjust = 0.5, 
    #   vjust = 0.5,
      angle = 45, hjust = 1
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