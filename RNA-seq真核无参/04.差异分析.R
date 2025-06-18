source("src/run_deseq.R")

run_deseq(
    fdata = "salmon_out/salmon_out.gene.counts.matrix",
    pdata = NULL,
    design_cols = "Group",
    contrast_list = c("Group", "B", "A"),
    pAdjustMethod = "fdr", # c("BH", "fdr", "none")
    pval_cutoff = 0.05,
    log2fc_cutoff = 1,
    smallestGroupSize = 3,
    outdir,
    threads = 8
)