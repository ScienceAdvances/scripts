using::using(using,ChIPseeker,GenomicFeatures,tidyverse,ggupset,clusterProfiler,forcats)
outdir <- 'Result/Enrich'
using::mkdir(outdir)
df <- data.table::fread('Result/ChIPseeker/WT_1_peakAnnotation.xls') %>% 
    dplyr::filter(distanceToTSS>= -3000|distanceToTSS<=3000)

colnames(df)
grep('TSS',df$annotation)
df$annotation %>% str_remove(' \\(.+\\)')
table(df$annotation %>% str_remove(' \\(.+\\)')) %>% as.data.frame()
geneid <- clusterProfiler::bitr(str_remove(df$geneId,'\\.\\d+'), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = 'org.Mm.eg.db')$ENTREZID

go_ea <- enrichGO(
    gene = geneid,
    #    universe      = names(geneList), # background
    OrgDb = 'org.Mm.eg.db',
    ont = "ALL", # BP", "MF", and "CC
    pAdjustMethod = "fdr",
    pvalueCutoff =1,
    qvalueCutoff = 1
)
GO <- as.data.frame(go_ea) %>% mutate() %>% group_by(ONTOLOGY) %>% slice_min(order_by=p.adjust,n=10)  %>% 
     dplyr::mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
head(GO,1)
colnames(GO)

p <- ggplot(GO) + 
  ggplot2::geom_bar(aes(FoldEnrichment,forcats::fct_reorder(Description,ONTOLOGY),fill = ONTOLOGY),stat="identity") + #将Term按照Ontology排序
  theme_bw() + 
  theme(text = element_text(family='serif',size =6,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black')) + 
#   geom_text(data = dplyr::filter(GO,qvalue<0.05) %>% dplyr::select(Count,Description),aes(x=Count+3,y=as.factor(Description),label='')) + #对于qvalue<0.05的term添加*
  labs(x='Fold Enrichment',y='') + 
  scale_fill_manual(values=c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD"))

using::gs(p,outdir = outdir,name='GO_ALL',w=6,h=9)

k <- clusterProfiler::enrichKEGG(
    gene = geneid,
    keyType = "kegg",
    organism = "mmu",
    pAdjustMethod = "fdr",
    pvalueCutoff = 1,
    qvalueCutoff = 1
)
# barplot
kegg <- as.data.frame(k) %>% slice_min(order_by=p.adjust,n=20,with_ties=FALSE)  %>% 
     dplyr::mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
kegg$Description %<>% stringr::str_remove(pattern = " - .+ \\(.+\\)")
colnames(kegg)

p <- ggplot(kegg) + 
  ggplot2::geom_bar(aes(FoldEnrichment,forcats::fct_reorder(Description,p.adjust),fill = p.adjust),stat="identity") + #将Term按照Ontology排序
  theme_bw() + 
  theme(text = element_text(family='serif',size =6,face = 'bold'),axis.text = element_text(family='serif',size =10,face = 'bold',colour = 'black')) + 
#   geom_text(data = dplyr::filter(GO,qvalue<0.05) %>% dplyr::select(Count,Description),aes(x=Count+3,y=as.factor(Description),label='')) + #对于qvalue<0.05的term添加*
  labs(x='Fold Enrichment',y='') + 
  ggplot2::scale_fill_gradient(high="#79B494",low="#D67E56")

using::gs(p,outdir = outdir,name='kegg',w=9,h=6)
