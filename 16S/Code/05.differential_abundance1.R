source('Code/src/RConfig.R')
outdir="Result/05.differential_abundance"
using::mkdir(outdir)
META=Sys.getenv("META")
# 创建microeco包可识别的整合对象
otu <- data.table::fread("Result/02.asv_taxonomy/ASV.xls", data.table = FALSE) %>%
  tibble::column_to_rownames("ASVID")
tax <- data.table::fread("Result/02.asv_taxonomy/taxonomy.xls", data.table = FALSE) %>%
  dplyr::select(-Confidence) %>%
  tibble::column_to_rownames("ASVID")

metadata <- data.table::fread(META, data.table = FALSE)
tree <- ape::read.tree("Result/02.asv_taxonomy/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

mt <- microtable$new(
  sample_table = metadata,
  otu_table = otu,
  tax_table = tax,
  phylo_tree = tree
)
mt %<>% tidy_taxonomy()
mt$tidy_dataset()
aaa=clone(mt)
# remove OTUs which are not assigned in the Kingdom “k__Archaea” or “k__Bacteria”.
mt$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”.
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# To make the OTU and sample information consistent across all files in the object
mt$tidy_dataset()
mt %<>% tidy_taxonomy()
# check the sequence numbers in each sample
mt$cal_abund(rel = TRUE)
## 根据最小测序深度进行数据抽平
set.seed(1314520)
mtable <- clone(mt) # 保留未抽平microtable数据对象。
mt_rarefied <- clone(mt)
mt_rarefied$rarefy_samples(sample.size = 10000)

conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::rename)
conflicts_prefer(lubridate::second)
conflicts_prefer(dplyr::desc)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::collapse)

group_order <- mt$sample_table$Group %>% unique()

taxa_level=c('all','Phylum','Famliy','Class','Genus','Species')
method=c("lefse", "rf", "metastat", "wilcox","DESeq2", "edgeR")
x='Phylum';y='lefse'

# test
for(x in taxa_level){
  mkdir(file.path(outdir,x))
  for(y in method){
  tmp <- clone(mt)
  if(!x =='all'){
    tmp$tax_table %<>% dplyr::filter(!!sym(x)!=glue::glue('{tolower(substr(x,1,1))}__'))
  }
  tmp$tidy_dataset()
  tmp$cal_abund()
  # tmp$sample_table$Group %<>% as.factor()
  t1 <- trans_diff$new(dataset = tmp, method =y,group = "Group", alpha = 1,remove_unknown=TRUE, taxa_level=x)
  t1$res_diff %>% as.data.table(keep.rownames='OTU') %>% mutate_if(is.numeric,round,digits=2) %>% 
    fwrite(file.path(outdir,glue::glue('{x}/{y}_table.xls')),sep='\t')
  tryCatch({
    p <- t1$plot_diff_bar(use_number = 1:10, width = 0.8, group_order = group_order,color_values = hue)
    using::gs(p, outdir = outdir, name = glue::glue('{x}/{y}_bar'), w = 9, h = 12)
    p <- t1$plot_diff_abund(plot_type = "barerrorbar", group_order = group_order, errorbar_addpoint = FALSE, errorbar_color_black = TRUE, plot_SE = TRUE,color_values = hue)
    using::gs(p, outdir = outdir, name = glue::glue('{x}/{y}_abundunce1'), w = 9, h = 12)
    p <- t1$plot_diff_abund(group_order = group_order, errorbar_addpoint = FALSE, errorbar_color_black = TRUE, plot_SE = TRUE,color_values = hue)
    using::gs(p, outdir = outdir, name = glue::glue('{x}/{y}_abundunce2'), w = 9, h = 12)
},
   error = function(e) e, finally = print(glue::glue('{x}/{y}')))
  }
}

rm(list=ls())

# remotes::install_github("lch14forever/microbiomeViz")
# ================ lefse ================
y='lefse'
t1 <- trans_diff$new(dataset = mt, method =y, group = "Group", alpha = 1, taxa_level='all', lefse_subgroup = NULL)
# The default is 4 which is showing the kingdom, Phylum, Class on the cladogram.
# taxa_level=list('Phylum'=5,'Class'=4,'Order'=3,'Famliy'=2,'Genus'=1,'Species'=0)
# for ( i in names(taxa_level)){

p <- t1$plot_diff_cladogram(
  color = hue200,
  filter_taxa = 0.0001, # 过滤掉相对丰度低于0.0001的分类单元。
  use_taxa_num = 200, # 根据丰度均值选择的展示分类单元。
  use_feature_num = 30, # 图中使用的特征数量,根据LDA得分筛选biomarker。
  clade_label_level = 5, # # 用字母标记标签的分类层级，5=Phylum，4=。
  only_select_show = FALSE, # 设置为TRUE，可以只展示select_show_labels选择的分类单元标签。
  group_order = group_order,
  node_size_offset = 1.5,
  annotation_shape = 21,# 处理图例不会一起变，仍是矩形。
  annotation_shape_size = 5,
  clade_label_size =3,
  alpha = 0.2
  )
using::gs(p, outdir = outdir, name = glue::glue('{y}_cladogram'), w = 16, h = 16)

# }
