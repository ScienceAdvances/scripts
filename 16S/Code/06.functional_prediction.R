source('Code/src/RConfig.R')
outdir <- 'Result/06.functional_prediction'
META=Sys.getenv("META")

using::using(ggpicrust2, tidyverse, ggprism, patchwork, readr, ggh4x, data.table, magrittr)
abundance <- data.table::fread("Result/06.functional_prediction/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", data.table = F) %>% column_to_rownames("function")
metadata <- data.table::fread(META, data.table = F)
# metadata$SampleID %<>% str_remove('_.+')
# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.
# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
d1 <- ggpicrust2::pathway_daa(abundance = abundance, p.adjust = "BH", metadata = metadata, group = "Group", daa_method = "ALDEx2", select = NULL, reference = NULL)

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
d2 <- d1 %>% dplyr::filter(method == "ALDEx2_Wilcoxon rank test")

# Annotate pathway results using KO to KEGG conversion
source('/home/data/wd/Vik/APP/ggpicrust2-main/R/pathway_annotation.R')
# a1=fread('/home/data/wd/Vik/DATAHUB/KEGG_orthology2pathway.txt')
# a2=fread('/home/data/wd/Vik/DATAHUB/pathway_annotation.txt')

# a3=merge(a1,a2,by.x='pathway_id',by.y='map_id') %>% distinct()
# fwrite(a3,'KEGG_Pathway.csv')
map <- fread('KEGG_Pathway.csv')
map$ko_id %<>% paste0('ko:',.)
d3 <- merge(map,d2,by.x='ko_id',by.y='feature') %>% distinct()
data.table::fwrite(d3,file.path(outdir,'KEGG.xls'),sep='\t')
head(d3,2)

c("pathway_name", "pathway_description", "pathway_class", "pathway_map")
d4 <- d3 %>% 
  dplyr::rename(feature=ko_id,pathway_map=pathway_id,pathway_class=PathwayLevel1,pathway_description=PathwayLevel2,pathway_name=PathwayLevel3)%>% 
  dplyr::arrange(p_values) %>% distinct(pathway_map,.keep_all=T) %>% distinct(feature,.keep_all=T)
# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(dplyr::summarise)
conflicted::conflicts_prefer(base::setdiff)
source('/home/data/wd/Vik/APP/ggpicrust2-main/R/pathway_errorbar.R')
kos <- d4 %>% dplyr::slice_min(as.numeric(p_values),n=20) %>% dplyr::pull(feature) %>% unique()
p <- pathway_errorbar(
  abundance = abundance,
  daa_results_df = as.data.frame(d4) %>% dplyr::filter(feature %in% kos),
  Group = metadata$Group,
  p_values_threshold = 0.05,
  order = "name",
  # select = kos,
  ko_to_kegg = F,
  p_value_bar = T,
  colors = hue,
  x_lab = "pathway_name"
) & theme_bw() & setheme()
using::gs(p, name = "KEGG", outdir = outdir, w = 20, h = 15)
