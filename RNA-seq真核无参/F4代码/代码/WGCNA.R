source("/home/data/wd/Vik/CODEHUB/RCode/config.R")
source('F4/src/wgcna_built_net.R')
source('F4/src/wgcna_picksoftpower.R')
source('F4/src/wgcna_qc.R')
source('F4/src/wgcna_select_modulegene.R')
outdir <- "F4/WGCNA"
suppressWarnings(dir.create(outdir, recursive = TRUE))
using::using(WGCNA, conflicted, impute, data.table, tidyverse)
fdata=fread('F2/数据/salmon.gene.TPM.not_cross_norm') %>% column_to_rownames('V1')
fdata <- log2(fdata+1)
rownames(fdata)  %>% unique %>% length
meta <- data.frame(
    Group=rep(c('CK5-6','CK7-8','CK9-10','NPA5-6', "NPA7-8", "NPA9-10"),each=3),
    row.names=colnames(fdata)
)
pdata <- model.matrix(~ Group - 1, data = meta) %>% as.data.frame()
colnames(pdata) %<>% str_remove('Group')

g1=fread('F2/差异分析/salmon.gene.counts.matrix.A1_vs_B1.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g2=fread('F2/差异分析/salmon.gene.counts.matrix.A2_vs_B2.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
g3=fread('F2/差异分析/salmon.gene.counts.matrix.A3_vs_B3.DESeq2.DE_results') %>% 
dplyr::filter(abs(log2FoldChange)>2,padj<0.01) %>% pull(1)
genes <- c(g1,g2,g3) %>% unique()
# venn_list <- list(g1,g2,g3)
# genes <- purrr::reduce(venn_list, base::intersect)

fdata <- fdata[genes,]
conflicted::conflict_prefer("filter", "dplyr")
qc_res <- wgcna_qc(
    fd = fdata %>% as.matrix(),
    pd = pdata,
    method = "average",
    # cutHeight = 166,
    cutMad = NULL,
    outdir = outdir,
    width = 26,
    height = 14
)
# Before qc: nSample = 18; nGenes = 34418
#  Flagging genes and samples with too many missing values...
#   ..step 1
# 没有缺失值，不进行删除

# After qc: nSample = 18; nGenes = 34418
softpower_res <- wgcna_picksoftpower(exp = qc_res$use_exp, outdir = outdir, net_type = "unsigned", Rsquared_cut = 0.9, width = 8, height = 6)
# The pick power: 7
conflicts_prefer(base::intersect)
wgcna_net <- wgcna_built_net(
    exp = qc_res$use_exp,
    cor_type="bicor",
    mergeCutHeight = 0.15,
    minModuleSize = 300,
    deepSplit=2,
    pheno = qc_res$use_pd,
    power = softpower_res$power,
    outdir = outdir,
    height = 9,
    width = 6
)
saveRDS(wgcna_net,file=str_glue("{outdir}/WGCNA.rds"))
wgcna_net <- readRDS(str_glue("{outdir}/WGCNA.rds"))
net=wgcna_net$net
write.table(wgcna_net$merged_infor, file = str_glue("{outdir}/merged_infor.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)

# 根据选定参数筛选WGCNA网络模块基因（核心基因），生成Figure_{select_pheno}_{select_module}_cor.pdf
merged_infor <- fread(str_glue("{outdir}/merged_infor.txt"), data.table = F)
conflict_prefer("filter", "dplyr")
# 筛选标准p<0.05且 |R|>0.3
# tan blue greenyellow midnightblue turquoise brown yellow
p <- fread('F4/WGCNA/SupplementaryTable_modTraitP.txt') %>% column_to_rownames('V1')
r <- fread('F4/WGCNA/SupplementaryTable_modTraitCor.txt') %>% column_to_rownames('V1')
colnames(merged_infor)
head(merged_infor)
conflict_prefer("filter", "dplyr")
merged_infor %>% head(1)
p
r
y='MEgreenyellow'
using::mkdir(file.path(outdir,'HubGene'))
for(y in setdiff(rownames(p),'MEgrey')){
    if(any(p[y,]<0.05 & abs(r[y,])>0.3)){
        color <- str_remove(y,'^ME')
        modulegene_res <- wgcna_select_modulegene(
            merged_infor = merged_infor,
            pheno_module_list = list(
                module = color,
                pheno = rep('CK5-6', 1)
            ),
            outdir = file.path(outdir,'HubGene')
        )
        hub_gene <- purrr::flatten_chr(modulegene_res) %>% unique()
        fwrite(data.table(Feature = hub_gene), file.path(outdir,'HubGene', str_glue("{color}.csv")))
    }
}


using::mkdir(file.path(outdir,'HubGene'))
for(x in colnames(pdata)){
    for(y in setdiff(rownames(p),'MEgrey')){
        if(p[y,x]<0.05 && abs(r[y,x])>0.3){
            color <- str_remove(y,'^ME')
            modulegene_res <- wgcna_select_modulegene(
                merged_infor = merged_infor,
                pheno_module_list = list(
                    module = color,
                    pheno = rep(x, 1)
                ),
                outdir = file.path(outdir,'HubGene')
            )
            hub_gene <- purrr::flatten_chr(modulegene_res) %>% unique()
            fwrite(data.table(Feature = hub_gene), file.path(outdir,'HubGene', str_glue("{x}_{color}.csv")))
        }
    }
}
