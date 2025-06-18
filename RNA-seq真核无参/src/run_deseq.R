#' @TODO 计算样本间差异--DESeq2 (测序数据)
#' @param countData 输入的原始read counts矩阵。行为基因，列为样本。
# countData可从以下工具获取
# package::function                     framework       output                DESeq2_input_function
# GenomicAlignments::summarizeOverlaps  R/Bioconductor  SummarizedExperiment  DESeqDataSet
# Rsubread::featureCounts               R/Bioconductor  matrix                DESeqDataSetFromMatrix
# tximport::tximport                    R/Bioconductor  list of matrices      DESeqDataSetFromTximport
# HTSeq::htseq-count                    Python          files                 DESeqDataSetFromHTSeq
# > head(countData)
#                    control1 control2 treat1 treat2
# ENSMUSG00000000001     1648     2306   2941   2780
# ENSMUSG00000000003        0        0      0      0
# ENSMUSG00000000028      835      950   1366   1051
# ENSMUSG00000000031       65       83     52     53
# ENSMUSG00000000037       70       53     94     66
# ENSMUSG00000000049        0        3      4      5
#' @param colData 输入的样本注释矩阵。colData的行名与countData的列名一致。
# because this table supplies metadata/information about the columns of the countData matrix.
# Notice that the first column of colData must match the column names of countData (except the first, which is the gene ID column).
# > head(colData)
#          condition  condition2
# control1   control  ttt
# control2   control  sss
# treat1       treat  sss
# treat2       treat  sss
#' @param design_cols 字符型向量，用于设定基于colData中的哪些列进行差异表达计算，可以是一个条件，也可以是多个
# design_cols="condition", design_formula= "~condition"
# design_cols=c("condition", "condition2"), design_formula= "~condition+condition2"
#' @param contrast 设定在计算差异表达的时候如何进行比较
# 当design_cols只有一个条件时，contrast是长度为3的字符型向量，分别指定
# the name of a factor in the design formula,
# the name of the numerator level for the fold change,
# and the name of the denominator level for the fold change (simplest case)
# 例如：
# contrast=c("condition", "A", "B") 结果中的foldchange = A/B
# contrast=c("condition", "B", "A") 结果中的foldchange = B/A
#' @param pval_cutoff padj阈值
#' @param log2fc_cutoff log2-transformed fold change阈值
#' @param return_df 是否以data.frame返回
#' @returnType
# return_df=TRUE  data.frame
# return_df=FALSE DESeqResults对象，类似于Tibbles，可以使用tidyverse包中的命令进行操作
#' @return 差异表达矩阵，如果没有则NA
#
#' @author Alex
#
run_deseq <- function(
    fdata = NULL,
    pdata = NULL,
    design_cols = "Group",
    contrast_list = c("Group", "B", "A"),
    pAdjustMethod = "fdr", # c("BH", "fdr", "none")
    pval_cutoff = 0.05,
    log2fc_cutoff = 1,
    smallestGroupSize = 3,
    outdir,
    threads = 8) {
    using::using(BiocParallel, DESeq2, tidyverse, data.table)
    # 构建DESeq2对象所需的design formula
    design_formula <- stats::as.formula(paste("~", paste(design_cols, collapse = "+")))
    # 构建DESeq2对象
    dds_obj <- DESeq2::DESeqDataSetFromMatrix(countData = fdata, colData = pdata, design = design_formula)
    # 基于DESeq2计算差异
    dds <- DESeq2::DESeq(dds_obj, parallel = TRUE, BPPARAM = BiocParallel::bpparam())
    keep <- rowSums(BiocGenerics::counts(dds) >= 10) >= smallestGroupSize
    dds <- dds[keep, ]
    # 提取结果
    lapply(contrast_list, function(x) {
        res <- DESeq2::results(dds, contrast = x, pAdjustMethod = pAdjustMethod, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = threads))
        data.table::as.data.table(as.data.frame(res), keep.rownames = "Feature") %>% data.table::fwrite(file.path(outdir, paste0(x, collapse = "_")))
    })
}
