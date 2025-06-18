source('/home/data/wd/Vik/CODEHUB/RCode/config.R')
using::using(data.table,tidyverse)
hue_day <- c('#F6E0D2','#DFA398','#9C6755','#659794','#EA967C','#F5C98E','#D65B5A','#586085')

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)

venn_list <- list(g1,g2,g3)
venn_gene <- purrr::reduce(venn_list, base::intersect)
names(venn_list) <-  c('CK5-6vsNPA5-6', "CK7-8vsNPA7-8", "CK9-10vsNPA9-10")
fwrite(data.table(GeneID = venn_gene), str_glue("F2/A交集基因.csv"))
# install.packages('ggvenn')
p <- ggvenn::ggvenn(venn_list,
    show_percentage = F, text_size = 8,
    fill_color = hue_day, fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
using::gs(p,outdir='F2',name='B')
