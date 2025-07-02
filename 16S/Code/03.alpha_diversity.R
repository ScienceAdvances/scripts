source("Code/src/RConfig.R")
using::using(picante)
outdir <- "Result/03.alpha_diversity"
using::mkdir(outdir)
hue <- paletteer::paletteer_c("scico::berlin", n = 150)
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
mpse %<>% MicrobiotaProcess::mp_cal_alpha(.abundance=RareAbundance)

f1 <- mpse %>% 
      MicrobiotaProcess::mp_plot_alpha(
        .group=Group, 
        .alpha=c('Observe', 'Chao1', 'ACE', 'Shannon', 'Simpson', 'Pielou')
      ) +
      scale_fill_manual(values=hue200, guide="none") +
      scale_color_manual(values=hue200, guide="none")+
        theme_bw() + setheme()
f2 <- mpse %>%
      MicrobiotaProcess::mp_plot_alpha(
        .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
      ) +
      scale_fill_manual(values=hue200, guide="none") +
      scale_color_manual(values=hue200, guide="none")+
      theme_bw() +
      setheme(45)
using::gs(f1 / f2 , outdir = outdir, name = 'AlphaDiversity', w = 28, h = 16)


# Faith’s index of PD遗传多样性指数-需要进化树才能计算。
RareAbundance=MicrobiotaProcess::mp_extract_assays(mpse, .abundance='RareAbundance', byRow = TRUE)

pd_result <- picante::pd(t(RareAbundance), tree, include.root=TRUE)
d <- merge(pd_result,metadata,by=0)

p <- ggpubr::ggboxplot(d,x='Group',y='PD',palette=hue200,fill='Group',add='jitter')+
  ggpubr::stat_compare_means()+setheme()
using::gs(p, outdir = outdir, name = 'FaithPD_index', w = 7, h = 6)

mpse %>% MicrobiotaProcess::mp_extract_sample()  %>% 
    dplyr::select(Sample,Observe, Chao1, ACE, Shannon, Simpson, Pielou) %>% 
    merge(pd_result,by.x='Sample',by.y=0) %>% 
    fwrite(file.path(outdir,'AlphaDiversity.xls'), sep='\t')