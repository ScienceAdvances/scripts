#' @TODO WGCNA样本和基因过滤
#' @title ### WGCNA样本和基因过滤
#' @description
#' @param fdata 表达谱，基因在行，样本在列
#' @param pdata 临床信息表，包含sample列，其它列为分析用的性状列，要求必须为数值型
#' @param method 样本聚类算法,默认average
#' @param cutHeight 需要剪切的枝长，默认为NULL，根据具体数值对样本进行筛选
#' @param cutMad 数值向量，mad中位数绝对偏差，默认为0.01；
#' 给定的cutMad和75%mad取最大值max_mad，保留中位数绝对偏差大于max_mad的基因；
#' 如果过滤以后基因数目过多，可以根据输出的mad信息调整cutMad值，top25保留的基因数目最少；
#' 如果为NULL，则不进行基因过滤
#' @param outdir 字符向量，输出路径
#' @param width 数值向量，生成图片的宽度
#' @param height 数值向量，生成图片的高度
#' @return
#' @examples qc_res <- wgcna_qc(
#'     fdata = train_data$tumor_fd,
#'     pdata = train_data$data_clinical %>% select(sample, Age) %>% mutate(Age = as.numeric(as.character(plyr::revalue(Age, c("Age<60" = 0, "Age>=60" = 1))))),
#'     method = "average", cutHeight = 120, cutMad = 0.01, outdir = outdir, width = 9, height = 6
#' )
#' @author *CY*
#'
wgcna_qc <- function(
    fdata = NULL,
    pdata = NULL,
    method = "average",
    cutHeight = NULL,
    cutMad = 0.01,
    outdir = base::getwd(),
    width = 9,
    height = 6) {
    # 输出路径
    using(tidyverse,WGCNA,fastcluster,data.table)
    mkdir(outdir)

    # 样本取交集
    common_sam <- base::intersect(colnames(fdata), rownames(pdata))
    if (length(common_sam) > 0) {
        fdata <- fdata[, common_sam, drop = FALSE]
        use_pd <- pdata[common_sam, , drop = FALSE]
        message(paste0("Before qc: nSample = ", ncol(fdata), "; nGenes = ", nrow(fdata)))

        # 如果mad不为空，根据mad中位数绝对偏差过滤基因；如果mad为空，则不过滤
        if (!is.null(cutMad)) {
            mad_res <- apply(fdata, 1, mad)
            max_mad <- max(quantile(mad_res, probs = seq(0, 1, 0.25))[2], cutMad) # 给定的cutMad和75%mad取最大值
            message(paste0("top 25% mad is ", sprintf("%.2f", quantile(mad_res, probs = seq(0, 1, 0.25))[4])))
            message(paste0("top 50% mad is ", sprintf("%.2f", quantile(mad_res, probs = seq(0, 1, 0.25))[3])))
            message(paste0("top 75% mad is ", sprintf("%.2f", quantile(mad_res, probs = seq(0, 1, 0.25))[2]), "; input cutMad is ", cutMad, "; max_mad is ", sprintf("%.2f", max_mad)))
            data_exp <- fdata[which(mad_res > max_mad), ]
        } else {
            data_exp <- fdata
        }

        # 转换为样品在行，基因在列的矩阵
        data_exp <- as.data.frame(t(data_exp))

        # 检测缺失值
        # 如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步。如果为FALSE，则需要进行删除缺失值。
        gsg <- goodSamplesGenes(data_exp, verbose = 3)
        if (!gsg$allOK) {
            # Optionally, print the gene and sample names that were removed:
            if (sum(!gsg$goodGenes) > 0) {
                printFlush(paste(
                    "Removing genes:",
                    paste(names(data_exp)[!gsg$goodGenes], collapse = ",")
                ))
            }
            if (sum(!gsg$goodSamples) > 0) {
                printFlush(paste(
                    "Removing samples:",
                    paste(rownames(data_exp)[!gsg$goodSamples], collapse = ",")
                ))
            }
            # Remove the offending genes and samples from the data:
            data_exp <- data_exp[gsg$goodSamples, gsg$goodGenes]
        } else {
            message("没有缺失值，不进行删除")
        }

        # 绘制样本聚类图
        # 根据cutHeight剔除离群样本
        sampleTree <- fastcluster::hclust(dist(data_exp), method = method)
        if (is.null(cutHeight)) {
            pdf(file = str_glue("{outdir}/Figure_SampleClustering_raw.pdf"), width = width, height = height)
            par(cex = 0.8)
            par(mar = c(3, 15, 10, 1))
            par(mgp = c(8, 1, 0))
            plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
            dev.off()
        } else {
            pdf(file = str_glue("{outdir}/Figure_SampleClustering_raw.pdf"), width = width, height = height)
            par(cex = 0.5)
            par(mar = c(3, 15, 10, 1))
            par(mgp = c(8, 1, 0))
            plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 5, cex.axis = 5, cex.main = 5)
            # 删除离群样本
            abline(h = cutHeight, col = "red") # 划定需要剪切的枝长
            dev.off()
            clust <- WGCNA::cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 1)
            keep_clust <- table(clust) %>%
                as.data.frame() %>%
                dplyr::slice_max(Freq,n=1) %>%
                dplyr::pull(clust) %>%
                as.character()
            data_exp <- data_exp[clust %in% keep_clust, ]
        }
        data.table::fwrite(as.data.table(data_exp, keep.rownames = "Sample"), file = str_glue("{outdir}/Table_WGCNA_fd.csv.gz"))
        use_pd[rownames(data_exp), , drop = FALSE] %>% 
            data.table::as.data.table(keep.rownames = 'Sample') %>% 
            data.table::fwrite(file = str_glue("{outdir}/Table_WGCNA_pd.csv"))
        message(paste0("\nAfter qc: nSample = ", nrow(data_exp), "; nGenes = ", ncol(data_exp)))
        return(list(use_exp = data_exp, use_pd = use_pd))
    } else {
        message("表达谱和临床信息的样本没有交集")
        return(NULL)
    }
}
