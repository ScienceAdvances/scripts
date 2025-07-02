source('Code/src/RConfig.R')
outdir <- "Result/07.statistics"
using::mkdir(outdir)
META=Sys.getenv("META")
using::using(igraph, ggraph, tidygraph, ggplot2, dplyr, magrittr)

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

ps_phylum <- phyloseq::tax_glom(ps, 'Genus')  # 聚合到门水平
otu <- phyloseq::otu_table(ps_phylum)
tax <- phyloseq::tax_table(ps_phylum)

d <- merge(tax[,'Genus',drop=F],otu,by=0)
d2 <- d %>% dplyr::filter(Genus!='g__') %>% 
    dplyr::mutate(Genus=str_remove(Genus,'g__')) %>%
    dplyr::select(-Row.names)

d2$Ref=apply(d2[,-1],1,sum)
e=d2 %>% dplyr::arrange(desc(Ref)) %>% 
    dplyr::distinct(Genus,.keep_all = T) %>%
    dplyr::slice_head(n=70) %>% 
    dplyr::select(-Ref)

fwrite(e,file.path(outdir,'Genus.xls'),sep='\t')


scriptdir=/home/data/wd/Vik/APP/SparCC3-master
# correlation calculation
python $scriptdir/SparCC.py \
    tmp/Genus.xls \
    --iter 20 \
    --algo SparCC \
    --thershold 0.1 \
    --cor_file tmp/cor_mat_sparcc.out

# pseudo p-value calculation
python $scriptdir/MakeBootstraps.py \
    tmp/Genus.xls \
    -n 5 \
    --path tmp/ \
    --template permutation_#.txt

python $scriptdir/SparCC.py tmp/permutation_0.txt --iter 5 --cor_file tmp/perm_cor_spp_0.txt 
python $scriptdir/SparCC.py tmp/permutation_1.txt --iter 5 --cor_file tmp/perm_cor_spp_1.txt  
python $scriptdir/SparCC.py tmp/permutation_2.txt --iter 5 --cor_file tmp/perm_cor_spp_2.txt 
python $scriptdir/SparCC.py tmp/permutation_3.txt --iter 5 --cor_file tmp/perm_cor_spp_3.txt 
python $scriptdir/SparCC.py tmp/permutation_4.txt --iter 5 --cor_file tmp/perm_cor_spp_4.txt

python $scriptdir/PseudoPvals.py \
    Result/07.statistics/cor_mat_sparcc.out \
    Result/07.statistics/perm_cor_spp_#.txt 5 \
    -o Result/07.statistics/pvals.xls \
    -t two_sided 


cor_sparcc <- read.delim(file.path(outdir,"cor_mat_sparcc.out"), row.names = 1, sep = '\t', check.names = FALSE)
#P值矩阵
pvals <- read.delim(file.path(outdir,"pvals.xls"), row.names = 1, sep = '\t', check.names = FALSE)

# 1. 数据准备 --------------------------------------------------------------
# 假设 SparCC 结果存储为三个文件：相关系数矩阵、p值矩阵、相对丰度表
cor_matrix <- as.matrix(cor_sparcc)
pval_matrix <- as.matrix(pvals)
abundance_df <- read.csv(file.path(outdir,"Genus.xls"),sep='\t', row.names = 1)

# 2. 数据筛选与转换 --------------------------------------------------------
# 将矩阵转换为边列表并筛选显著相关对
edges <- expand.grid(Genus1 = rownames(cor_matrix), Genus2 = colnames(cor_matrix)) %>%
  dplyr::filter(as.numeric(Genus1) < as.numeric(Genus2)) %>%  # 去除重复对和自相关
  dplyr::mutate(
    cor = cor_matrix[cbind(Genus1, Genus2)],
    pval = pval_matrix[cbind(Genus1, Genus2)]
  ) %>%
  dplyr::filter(abs(cor) > 0.3, pval < 0.05)  # 按图中阈值筛选
head(edges,10)
glimpse(edges)
# 3. 节点属性准备 ----------------------------------------------------------
nodes <- data.frame(
  Genus = rownames(abundance_df),
  Abundance = rowMeans(abundance_df)  # 计算平均相对丰度
) %>%
  dplyr::filter(Genus %in% unique(c(edges$Genus1, edges$Genus2)))  # 仅保留网络中的节点

# 4. 构建网络图对象 --------------------------------------------------------
head(edges,2)
edges$Genus1 %<>% as.character()
graph <- as_tbl_graph(edges) %>%
  activate(nodes) %>%
  dplyr::left_join(nodes, by = c("name" = "Genus")) %>%
  dplyr::mutate(
    size = sqrt(Abundance) * 30,  # 节点大小与丰度平方根成比例
    label = if_else(Abundance > 0.01, name, "")  # 仅标注丰度>1%的属
  )

# 5. 可视化参数设置 --------------------------------------------------------
set.seed(123)  # 固定布局
layout_type <- create_layout(graph, layout = "fr")  # Fruchterman-Reingold布局

# 6. 绘制高级网络图 --------------------------------------------------------
p <- ggraph(layout_type) +
  geom_edge_link(
    aes(
      edge_width = abs(cor),  # 边宽=相关系数绝对值
      edge_color = ifelse(cor > 0, hue4[2],hue4[1])  # 红正绿负
    ),
    alpha = 0.8, 
    arrow = arrow(length = unit(2, "mm")), 
    end_cap = circle(3, "mm")
  ) +
  geom_node_point(
    aes(size = size),
    color = "steelblue",
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_node_text(
    aes(label = label),
    size = 3.5,
    repel = TRUE,  # 防标签重叠
    color = "black",
    bg.color = "white",
    bg.r = 0.15
  ) +
  scale_edge_width_continuous(
    name = "|Correlation|",
    range = c(0.5, 3),  # 边宽范围
    breaks = c(0.3, 0.6, 0.9)  # 图例刻度
  ) +
  scale_edge_color_manual(
    name = "Direction",
    values = c(hue4[2], hue4[1]),
    labels = c("Negative", "Positive")
  ) +
  guides(
    edge_color = guide_legend(override.aes = list(edge_width = 3)),
    edge_width = guide_legend(override.aes = list(edge_color = "grey50"))
  ) +
  theme_graph(base_family = "sans") +
  labs(
    title = "Genus-level Co-occurrence Network (SparCC)",
    subtitle = "|Correlation| > 0.3 & p < 0.05\nNode size ~ Relative Abundance",
    caption = "SparCC algorithm with FDR correction"
  )

# 7. 保存高清图片 ---------------------------------------------------------
using::gs(p,format=c('pdf','png'),name='network_plot',outdir = outdir,w=12,h=10)
