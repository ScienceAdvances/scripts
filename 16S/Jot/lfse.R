source("/home/data/wd/Vik/X155/Code/src/RConfig.R")
outdir <- "Result/05.Model"
using::mkdir(outdir)
library(MicrobiotaProcess)
color2 <- hue4[2:1]

otu <- data.table::fread("Result/dada2/ASV.xls", data.table = FALSE) %>%
  tibble::column_to_rownames("ASVID")
tax <- data.table::fread("Result/dada2/taxonomy.xls", data.table = FALSE) %>%
  dplyr::select(-Confidence) %>%
  tibble::column_to_rownames("ASVID")

metadata <- data.table::fread("Result/metadata.tsv", data.table = FALSE)
tree <- ape::read.tree("Result/dada2/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(otu),taxa_are_rows=TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  sample_data(metadata),
  phy_tree(tree))
mpse <- MicrobiotaProcess::as.MPSE(ps)

mpse %<>% MicrobiotaProcess::mp_rrarefy() %>% 
MicrobiotaProcess::mp_cal_abundance(.abundance = RareAbundance) %>% 
mpse %<>% MicrobiotaProcess::mp_cal_abundance(.abundance = RareAbundance,.group=Group)
RareAbundance=mp_extract_assays(mpse, .abundance='RareAbundance', byRow = TRUE)

mpse %<>%
    MicrobiotaProcess::mp_diff_analysis(
        .abundance = RelRareAbundanceBySample,
        .group = Group,
        first.test.alpha = 0.05,
        filter.p="pvalue"
    )

taxa.tree <- mpse %>% 
  mp_extract_tree(type="taxatree")
taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, pvalue, fdr) %>% dplyr::filter(!is.na(fdr)) %>% dplyr::arrange(pvalue)

p1 <- ggtree::ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  )

p2 <- p1 +
 ggtree::geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  )

library(ggtreeExtra)
library(ggstar)
p3 <- p2 +
  ggnewscale::new_scale("fill") +
  ggtreeExtra::geom_fruit(
    data = MicrobiotaProcess::td_unnest(RareAbundanceBySample),
    geom = ggstar::geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Group, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Group,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=color2)

p4 <- p3 + geom_tiplab(size=2, offset=7.2)
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_time,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_time)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_time,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02"))

p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)
p6