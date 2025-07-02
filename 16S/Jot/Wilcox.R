# 加载包
library(dplyr)
library(tidyr)
library(broom)
library(vegan)

# 假设输入数据格式：
# - otu_table: 行为样本，列为物种（如Genus），数值为丰度（绝对或相对）
# - metadata: 包含分组信息（如Group列，取值为"HC"和"UC"）
otu <- data.table::fread("Result/02.稀释曲线/otu.xls",data.table=FALSE) %>% tibble::column_to_rownames("ASV") %>% as.matrix()
tax <- data.table::fread("Result/02.稀释曲线/tax.xls",data.table=FALSE) %>% tibble::column_to_rownames("ASV") %>% as.matrix()
metadata <- data.table::fread("Result/metadata.tsv",data.table=FALSE)
tree <- ape::read.tree("Result/01.Report/rooted_tree.nwk")
rownames(metadata) <- metadata$SampleID

ps <- phyloseq(
  phyloseq::otu_table(otu,taxa_are_rows=TRUE),
  phyloseq::tax_table(tax),
  sample_data(metadata),
  phy_tree(tree))
ps = prune_taxa(taxa_sums(ps) > 1, ps)


# 聚合到种水平（需确保 tax_table 中存在 Species 列）
physeq_species_agg <- tax_glom(ps, taxrank = "Genus")

# 提取丰度表
species_abundance <- otu_table(physeq_species_agg) %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU_ID")

# 提取分类表
species_taxonomy <- phyloseq::tax_table(physeq_species_agg) %>%
  as.data.frame() %>%
  select(Phylum, Genus)  # 选择需要的分类层级

# 合并丰度与分类信息
species_table <- species_abundance %>%
  left_join(rownames_to_column(species_taxonomy,'OTU_ID'), by = "OTU_ID" ) %>%
  relocate(Phylum, Genus, .after = OTU_ID)

species_table[1:3,1:2]
# 保存为 CSV
write.csv(species_table, "species_level_otu_table.csv", row.names = FALSE)

# 定义函数：计算每个物种的统计量及p值
calculate_metastats <- function(otu, group) {
  # 按分组计算均值、方差、标准误
  group_stats <- data.frame(
    Mean_HC = mean(otu[group == "HC"]),
    Variance_HC = var(otu[group == "HC"]),
    Std.Err_HC = sd(otu[group == "HC"]) / sqrt(sum(group == "HC")),
    Mean_UC = mean(otu[group == "UC"]),
    Variance_UC = var(otu[group == "UC"]),
    Std.Err_UC = sd(otu[group == "UC"]) / sqrt(sum(group == "UC"))
  )
  
  # 非参数检验（Mann-Whitney U检验，替代Metastats的非参数t检验）
  p_value <- wilcox.test(otu ~ group, exact = FALSE)$p.value
  
  list(stats = group_stats, p_value = p_value)
}

x=t(otu_table)[,1]
b=species_table[,-c(1,2)]
b$Ref=apply(b[,-1],1,sum)
e=b %>% dplyr::arrange(desc(Ref)) %>% dplyr::distinct(Genus,.keep_all = T) %>% dplyr::select(-Ref)

# 遍历所有物种进行计算
results <- lapply(e %>% column_to_rownames('Genus') %>% t() %>% as.data.frame(), function(x) {
  res <- calculate_metastats(x, metadata$Group)
  data.frame(res$stats, p_value = res$p_value)
})

# 合并结果
result_table <- bind_rows(results, .id = "Genus") %>%
  dplyr::mutate(
    p_value = round(p_value, 4),
    q_value = p.adjust(p_value, method = "fdr")  # FDR校正[1](@ref)
  ) %>%
  dplyr::arrange(p_value)

result_table$Genus %<>% str_remove('g__')
result_table %>% mutate_if(is.numeric,round,digits=2) %>% 
fwrite('Table/Metastats.xls',sep='\t')
