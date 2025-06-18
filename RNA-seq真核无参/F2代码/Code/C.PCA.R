source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(data.table,tidyverse)
hue_day <- c('#F6E0D2','#DFA398','#9C6755','#659794','#EA967C','#F5C98E','#D65B5A','#586085')

d=fread('F2/数据/salmon.gene.TPM.not_cross_norm') %>% column_to_rownames('V1')
d=log2(d+1)
gr <- rep(c('CK5-6','CK7-8','CK9-10','NPA5-6', "NPA7-8", "NPA9-10"),each=3)

install.packages('factoextra')
using::using(FactoMineR,factoextra)

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
genes <- c(g1,g2,g3) %>% unique()


iris.pca = PCA(t(d[genes,]), scale.unit = TRUE, ncp = 2, graph = F)
p=factoextra::fviz_pca_ind(X = iris.pca, ## pca对象
  axes = 1:2, ## 展示的两个主成分
  geom = 'point', ## 展示individual的形式
  habillage = as.factor(gr), #分组 factor
  legend.title = 'Groups', ## 分组变量的title
  palette = hue_day, ## 颜色面板
  addEllipses = F,  ## 是否绘制椭圆
#   ellipse.level = 0.9, ## 椭圆的大小
  title = 'PCA Plot', ## 标题
  mean.point = F ## 不删除每个组的重心
  )+
  theme(plot.title = element_text(hjust = 0.5,size=18))
using::gs(p,outdir='F2',name='C',w=4,h=3)
