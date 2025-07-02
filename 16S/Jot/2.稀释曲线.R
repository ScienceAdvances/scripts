source("Code/src/RConfig.R")
outdir <- "Result/02.稀释曲线"
using::mkdir(outdir)

# 创建microeco包可识别的整合对象
otu <- data.table::fread("Result/01.Report/ASV.xls", data.table = FALSE) %>%
  tibble::column_to_rownames("ASVID") %>%
  as.data.frame()
tax <- data.table::fread("Result/01.Report/taxonomy.xls", data.table = FALSE) %>%
  dplyr::select(-Confidence) %>%
  tibble::column_to_rownames("ASVID") %>%
  as.data.frame()

metadata <- data.table::fread("Result/metadata.tsv", data.table = FALSE)
tree <- ape::read.tree("Result/01.Report/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

dataset <- microtable$new(
  sample_table = metadata,
  otu_table = otu,
  tax_table = tax,
  phylo_tree = tree
)
dataset %<>% tidy_taxonomy()

# remove OTUs which are not assigned in the Kingdom “k__Archaea” or “k__Bacteria”.
dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast”.
dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# To make the OTU and sample information consistent across all files in the object
dataset$tidy_dataset()
# check the sequence numbers in each sample
dataset$cal_abund()

# 1.2.3 数据抽平-mecodev包
## 稀释度度曲线
dataset$sample_sums() %>% sort() # 计算每个样品种包含序列数。# 升序排列
set.seed(1314520)
t1 <- trans_rarefy$new(
  dataset,
  alphadiv = "Shannon",
  depth = c(0, 10, 50, 500, seq(1000, max(dataset$sample_sums()), 1000))
)
t1$res_rarefy %>% head() # 查看稀释度采样结果


str(dataset, 1)
dataset$sample_table$Group %>% table()
## 绘制稀释度曲线图
p <- t1$plot_rarefy(
  color = "SampleID",
  color_values = hue4,
  show_point = F,
  add_fitting = F
) +
  theme_bw()+setheme()
using::gs(p, outdir = outdir, name = "rarefy-Sample", w = 12, h = 8)

p <- t1$plot_rarefy(
  color = "Group",
  color_values = hue4
) +
  theme_bw()+setheme()
using::gs(p, outdir = outdir, name = "rarefy-Group", w = 12, h = 8)
## 看图可以按照最低测序样本的测序深度进行抽平。

## 根据最小测序深度进行数据抽平
set.seed(1314520)
data <- clone(dataset) # 保留未抽平microtable数据对象。
data$rarefy_samples(sample.size = min(dataset$sample_sums()))
data$sample_sums() %>% range() # 计算每个样品种包含序列数。 # 最小值与最大值

otu_flattenning <- data$otu_table
tax_flattenning <- data$tax_table
fwrite(as.data.table(otu_flattenning, keep.rownames = "ASV"), "Result/02.稀释曲线/otu.xls", sep = "\t")
fwrite(as.data.table(tax_flattenning, keep.rownames = "ASV"), "Result/02.稀释曲线/tax.xls", sep = "\t")

# 1.2.4 指定分类水平顺序
data$sample_table$Group %<>% factor(., levels = unique(metadata$Group))

# 特征统计表
data.table(
  SampleID=colnames(otu_flattenning),
  Number_of_Features=apply(otu_flattenning,2,function(x){sum(x!=0)}),
  Number_of_Sequences=apply(otu_flattenning,2,sum)
  ) %>% fwrite('Result/01.Report/Table2.xls',sep='\t')

# 读取OTU表（行为样本，列为物种）
otu <- read.table("otu.txt", sep = "\t", row.names = 1, header = TRUE)
otu <- otu[, -ncol(otu)]  # 删除分类信息列
otu <- t(otu)  # 转置为vegan要求的格式（行为样本）

# # 直接绘制稀释曲线（默认Richness指数）
# pdf('X.pdf')
# vegan::rarecurve(otu, 
#   sample='Group',
#   step = 1000,  # 抽样步长（根据测序量调整）
#   xlab = "Sequencing Depth (Reads)", 
#   ylab = "Observed Species Richness",
#   label = F,
#   col=hue4
# )  # 避免样本标签重叠
# dev.off()