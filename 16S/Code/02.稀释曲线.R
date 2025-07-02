source("Code/src/RConfig.R")
library(MicrobiotaProcess)
outdir <- "Result/07.statistics"
using::mkdir(outdir)
hue <- c(hue4,paletteer::paletteer_c("scico::berlin", n = 150))
META=Sys.getenv("META")

otu <- data.table::fread("Result/02.asv_taxonomy/ASV.xls", data.table = FALSE) %>%
  tibble::column_to_rownames("ASVID")
tax <- data.table::fread("Result/02.asv_taxonomy/taxonomy.xls", data.table = FALSE) %>%
  dplyr::select(-Confidence) %>%
  tibble::column_to_rownames("ASVID")

metadata <- data.table::fread(META, data.table = FALSE)
tree <- ape::read.tree("Result/02.asv_taxonomy/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(otu),taxa_are_rows=TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  sample_data(metadata),
  phy_tree(tree))
mpse <- MicrobiotaProcess::as.MPSE(ps)

mpse %<>% MicrobiotaProcess::mp_rrarefy() %>% MicrobiotaProcess::mp_cal_abundance(.abundance = RareAbundance)
RareAbundance=mp_extract_assays(mpse, .abundance='RareAbundance', byRow = TRUE)

data.table::as.data.table(RareAbundance, keep.rownames = "ASVID") %>% 
  data.table::fwrite(file.path(outdir,'抽平后的丰度表.xls'), sep = "\t")

# 特征统计表
data.table(
  SampleID=colnames(RareAbundance),
  Number_of_Features=apply(RareAbundance,2,function(x){sum(x!=0)}),
  Number_of_Sequences=apply(RareAbundance,2,sum)
) %>% fwrite(file.path(outdir,'Feature_Stat.xls'),sep='\t')

relative_abundance <- function(table){
    tibble::as_tibble(table,rownames='ASV') %>% 
        tidyr::pivot_longer(cols=-ASV) %>% 
        dplyr::group_by(name) %>% 
        dplyr::mutate(total=sum(value)) %>% 
        dplyr::mutate(relative_abun=value/total) %>% 
        tidyr::pivot_wider(id_cols=c("ASV"),names_from = name,values_from = relative_abun) %>% 
        dplyr::mutate_if(is.numeric,round,digits=5)
}
# 相对丰度表
data.table::fread("Result/vsearch/ASV.xls", data.table = FALSE) %>%
    tibble::column_to_rownames("ASVID") %>% 
    relative_abundance() %>% 
    data.table::fwrite(file.path(outdir,'相对丰度表.xls'),sep='\t')

data.table::fread(file.path(outdir,'抽平后的丰度表.xls'), data.table = FALSE) %>%
    tibble::column_to_rownames("ASVID") %>% 
    relative_abundance() %>% 
    data.table::fwrite(file.path(outdir,'抽平后的相对丰度表.xls'), sep='\t')

# ============= 稀释曲线 ==============
mpse %<>% MicrobiotaProcess::mp_cal_rarecurve(.abundance = RareAbundance, chunks = 400)

alphas <- c('Chao1','ACE','Shannon')
p <- mpse %>% MicrobiotaProcess::mp_plot_rarecurve(
      .rare = RareAbundanceRarecurve, 
      .alpha = !!alphas,
      .group = 'Group', 
      plot.group = TRUE
  ) +
  scale_color_manual(values=hue) +
  scale_fill_manual(values=hue,guide="none")  +
  labs(x='Number of Sequences Sampled',y='Alpha Diversity Index') + theme_bw()+
  setheme() 
using::gs(p, outdir = outdir, name = 'AlphaDiversityIndex-Group-RareCurve', w = 18, h = 6)
p <- mpse %>% MicrobiotaProcess::mp_plot_rarecurve(.rare = RareAbundanceRarecurve, .alpha = !!alphas,nrow=1)+
  scale_color_manual(values=hue) +
  scale_fill_manual(values=hue,guide="none")  +
  labs(x='Number of Sequences Sampled',y='Alpha Diversity Index')+ 
  theme_bw() + setheme()
using::gs(p, outdir = outdir, name = 'AlphaDiversityIndex-Sample-RareCurve', w = 18, h = 6)

p <-mpse %>% MicrobiotaProcess::mp_plot_rarecurve(
      .rare = RareAbundanceRarecurve, 
      .alpha = Observe,
      .group = 'Group', 
      plot.group = TRUE
  ) +
  scale_color_manual(values=hue) +
  scale_fill_manual(values=hue,guide="none")  +
  labs(x='Number of Sequences Sampled',y='Number of Features') + theme_bw()+
  setheme() 
using::gs(p, outdir = outdir, name = 'RareCurve-Group', w = 9, h = 6)
p <- mpse %>% MicrobiotaProcess::mp_plot_rarecurve(.rare = RareAbundanceRarecurve, .alpha = Observe)+
  scale_color_manual(values=hue) +
  scale_fill_manual(values=hue,guide="none")  +
  labs(x='Number of Sequences Sampled',y='Number of Features') + 
  theme_bw() + setheme()
using::gs(p, outdir = outdir, name = 'RareCurve-Sample', w = 9, h = 6)

# ============= RankAbundanceCurves ==============
myrank <- function(x){
  res <- length(x)-rank(x, ties.method = "first")+1
  res[x==0] <- NA
  return(res)
}

# 生成rank矩阵
otutab_rank <- map_dfr(otu, myrank) %>% as.data.frame()
rownames(otutab_rank) <- rownames(otu)

# 重塑数据
otutab2 <- otu %>% rownames_to_column("FeatureID") %>% gather("Sample", "Count", -FeatureID)
otutab_rank2 <- otutab_rank %>% rownames_to_column("FeatureID") %>%
  gather("Sample", "Rank", -FeatureID)

# 合并otutab数据和ranks数据
data <- left_join(otutab2, otutab_rank2) %>% group_by(Sample) %>%
  dplyr::mutate(Abundance = Count/sum(Count)) %>% na.omit()
df <- merge(data, metadata , by.x = 'Sample', by.y = 0)

head(otutab2,2)
head(otutab_rank2,2)
# 绘制 Rank-Abundance，注意y轴数据进行对数变化

p <- ggplot(df, aes(x = Rank, y = Abundance, color = Group)) + 
  geom_line(aes(group=Sample)) +
  scale_y_log10() +
  scale_color_manual(values=hue) +
  labs(title="Rank-Abundance Curves") + ylab("% Abundance")+
  theme_bw() +
  setheme()
using::gs(p, outdir = outdir, name = 'RankAbundanceCurves', w = 7, h = 6)
