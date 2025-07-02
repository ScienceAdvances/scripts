# 安装picrust2
micromamba env create -f /home/data/wd/Vik/APP/picrust2-master/picrust2-env.yaml
micromamba clean --all
micromamba activate picrust2
cd /home/data/wd/Vik/APP/picrust2-master
pip install --editable .

chmod u+x -R /home/data/wd/Vik/APP/picrust2-master/picrust2
ll /home/data/wd/Vik/APP/picrust2-master/picrust2
cd /home/data/wd/Vik/X155
rm -rf Result/picrust2
/home/data/wd/Vik/APP/picrust2-master/picrust2/MinPath/MinPath12hmp.py -h
picrust2_pipeline.py \
--study_fasta Result/vsearch/RepresentSequence.fasta \
--input Result/vsearch/ASV.xls \
--output Result/picrust2 \
--min_align 0.8 \
--processes 16

# lapply(list('DESeq2','R.utils','edgeR','ALDEx2','lefser','Maaslin2','metagenomeSeq'),BiocManager::install)
# 安装ggpicrust2
source('Code/src/RConfig.R')
outdir <- 'Result/picrust2'
# remotes::install_github("cafferychen777/ggpicrust2")
using::using(ggpicrust2, tidyverse, ggprism, patchwork, readr, ggh4x, data.table, magrittr)

abundance <- data.table::fread("Result/picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", data.table = F) %>% column_to_rownames("function")
metadata <- data.table::fread("Result/metadata.tsv", data.table = F)

# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.
# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
d1 <- ggpicrust2::pathway_daa(abundance = abundance, p.adjust = "BH", metadata = metadata, group = "Group", daa_method = "DESeq2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
d2 <- d1 %>% dplyr::filter(method == "ALDEx2_Wilcoxon rank test")
d2 <- d1
d1 %>% head(2)
d1$method %>% unique()

# Annotate pathway results using KO to KEGG conversion
# source('/home/data/wd/Vik/APP/ggpicrust2-main/R/pathway_annotation.R')
# d3 <- pathway_annotation(pathway = "KO", daa_results_df = d2, ko_to_kegg = TRUE,p_values_threshold=0.05)
# data.table::fwrite(d3,file.path(outdir,'KEGG.xls'),sep='\t')
head(d3,2)
# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
source('/home/data/wd/Vik/APP/ggpicrust2-main/R/pathway_errorbar.R')
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarise)
p <- pathway_errorbar(
  abundance = abundance,
  daa_results_df = d3,
  Group = metadata$Group,
  p_values_threshold = 1,
  order = "group",
  select = d3 %>% dplyr::slice_min(p_values,n=20) %>% dplyr::pull(feature),
  ko_to_kegg = F,
  p_value_bar = T,
  colors = hue4[2:1],
  x_lab = "pathway_name"
) & theme_bw() & setheme()
using::gs(p, name = "KO", outdir = "Result/picrust2", w = 25, h = 18)

head(daa_annotated_sub_method_results_df)
daa_annotated_sub_method_results_df$pathway_class %>%
  table() %>%
  data.frame()
# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.


abundance <- data.table::fread("Result/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", data.table = F) %>% column_to_rownames("function")
metadata <- data.table::fread("Result/metadata.tsv", data.table = F)
# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
e1 <- pathway_daa(abundance = abundance, metadata = metadata, group = "Group", p.adjust = "BH", daa_method = "ALDEx2", select = NULL, reference = NULL)
head(e1,2)
# Filter results for ALDEx2_Kruskal-Wallace test method
e2 <- e1 %>% dplyr::filter(method == "ALDEx2_Wilcoxon rank test")
head(daa_annotated_sub_method_results_df,2)
daa_results_df$method %>% unique()
# Annotate pathway results without KO to KEGG conversion
e3 <- pathway_annotation(pathway = "EC", daa_results_df = e2, ko_to_kegg = FALSE)
e3$method %<>% as.character()
# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset


source('/home/data/wd/Vik/APP/ggpicrust2-main/R/pathway_errorbar.R')
p <- pathway_errorbar(
  abundance = abundance,
  daa_results_df = e3,
  Group = metadata$Group,
  p_values_threshold = 1,
  order = "group",
  select = e3 %>% dplyr::slice_min(p_values,n=20) %>% dplyr::pull(feature),
  ko_to_kegg = FALSE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
) & theme_bw()
using::gs(p, name = "EC", outdir = "Result/ggpicrust2", w = 9, h = 9)

head(daa_results_df,2)