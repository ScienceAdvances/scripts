source('Code/src/RConfig.R')
using::using(ggprism,MicrobiotaProcess)
library(MicrobiotaProcess)
outdir="Result/04.beta_diversity"
META=Sys.getenv("META")

otu <- data.table::fread("Result/02.asv_taxonomy/ASV.xls", data.table = FALSE) %>%
  tibble::column_to_rownames("ASVID")
tax <- data.table::fread("Result/02.asv_taxonomy/taxonomy.xls", data.table = FALSE) %>%
  dplyr::select(-Confidence) %>%
  tibble::column_to_rownames("ASVID")

metadata <- data.table::fread(META, data.table = FALSE)
tree <- ape::read.tree("Result/02.asv_taxonomy/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

mpse <- phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(otu),taxa_are_rows=TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  sample_data(metadata),
  phy_tree(tree)) %>% 
  MicrobiotaProcess::as.MPSE()

mpse %<>% MicrobiotaProcess::mp_rrarefy() %>% MicrobiotaProcess::mp_cal_abundance(.abundance = RareAbundance)
mpse %<>% mp_decostand(.abundance=Abundance, method = "hellinger", logbase = 2)
mpse %<>% MicrobiotaProcess::mp_cal_alpha(.abundance=RareAbundance)

# bray
mpse %<>% mp_cal_dist(.abundance = 'hellinger', distmethod = "bray") %>% 
    mp_adonis(.abundance='hellinger', .formula=~Group, distmethod="bray", permutations=9999, action="add") %>% 
    mp_cal_pca(.abundance='hellinger', distmethod= 'bray') %>% 
    mp_cal_pcoa(.abundance='hellinger', distmethod= 'bray') %>% 
    mp_cal_nmds(.abundance='hellinger', distmethod= 'bray')

betaplot(mpse,ordination='PCA',distance='bray',hue=hue200,w=16,h=7)
betaplot(mpse,ordination='PCoA',distance='bray',hue=hue200,w=16,h=7)
betaplot(mpse,ordination='NMDS',distance='bray',hue=hue200,w=16,h=7)

# wunifrac
for( i in c('unifrac','jaccard','wunifrac')){
  mpse %<>% mp_cal_dist(.abundance = 'RareAbundance', distmethod = i, action="add") %>% 
      mp_adonis(.abundance='RareAbundance', .formula=~Group, distmethod=i, permutations=9999, action="add") %>% 
      mp_cal_pca(.abundance='RareAbundance', distmethod= i) %>% 
      mp_cal_pcoa(.abundance='RareAbundance', distmethod= i) %>% 
      mp_cal_nmds(.abundance='RareAbundance', distmethod= i)
  betaplot(mpse,ordination='PCA',distance=i,hue=hue200,w=16,h=7)
  betaplot(mpse,ordination='PCoA',distance=i,hue=hue200,w=16,h=7)
  betaplot(mpse,ordination='NMDS',distance=i,hue=hue200,w=16,h=7)
}
walk(c("bray",'unifrac','jaccard','wunifrac'),function(x){
  mpse %>% mp_extract_dist(distmethod=x) %>% as.matrix() %>% as.data.table(keep.rownames='SampleID') %>% 
  fwrite(file.path(outdir,glue::glue("{x}_BetaDiversity.xls")),sep='\t')
})
RareAbundance=mp_extract_assays(mpse, .abundance='RareAbundance', byRow = TRUE)

# 应力曲线 应力函数值（<=0.2合理）
# 计算Bray-Curtis相异性矩阵
dist <- vegan::vegdist(RareAbundance %>% t(), method = "bray")
# 执行NMDS分析
swiss_mds <- vegan::metaMDS(dist)
#展示结果
summary(swiss_mds)
pdf(file.path(outdir,"stress.pdf"))
stressplot(swiss_mds, p.col = hue4[1], l.col = hue4[2])
dev.off()


mpse %<>% mp_adonis(.abundance=hellinger, .formula=~Group, distmethod="bray", permutations=9999, action="add")


# =================== anosim ===================
df_anosim <- vegan::anosim(RareAbundance %>% t(), grouping=factor(metadata$Group), permutations = 999,parallel = getOption("mc.cores"))#数据也可以是原始otu数据
#整理出作图数据
df1 <- data.frame(
  x=df_anosim$class.vec,
  y=df_anosim$dis.rank
)

#绘图
p <- ggplot(df1,aes(x=x,y=y))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=x), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  scale_fill_manual(values=hue200)+ #指定颜色
  ggtitle("Bray-Curtis Anosim")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 14,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+
  theme(legend.position = 'none')+
  labs(x = paste("R=",df_anosim$statistic,", ","p=", df_anosim$signif),
       y = "Rank of Distance (Bray_Curtis)")
using::gs(p, outdir = outdir, name = 'anosim', w = 7, h = 6)