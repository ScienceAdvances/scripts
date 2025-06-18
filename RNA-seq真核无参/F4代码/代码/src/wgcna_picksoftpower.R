#' @TODO 筛选用于构建wgcna网络的软阈值
#' @title ### 筛选用于构建wgcna网络的软阈值
#' @param exp 表达谱文件（wgcna_qc生成的表达谱）样本在行，基因在列
#' @param outdir 字符向量，输出路径
#' @param net_type 字符向量,构建网络的类型;默认unsigned，可选 "signed" 或 "signed hybrid"
#' @param Rsquared_cut 数值，默认0.85。筛选最优软阈值的临界值，最低标准是0.8
#' @param width 数值向量，生成图片的宽度
#' @param height 数值向量，生成图片的高度
#' @return *list*
#' @examples softpower_res <- pickSoftPower(exp = qc_res$use_exp, outdir = outdir, net_type = "unsigned", Rsquared_cut = 0.85, width = 8, height = 6)
#' @author *CY*
#'
wgcna_picksoftpower <- function(exp = NULL, outdir = NULL, net_type = "unsigned",
                                Rsquared_cut = 0.85, width = 8, height = 6) {
    mkdir(outdir)
    # 输出样本和基因数目
    nSamples <- nrow(exp)
    message(paste0("Input Gene: ", ncol(exp), ";Input Sample: ", nSamples))

    # 设定软阈值范围
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    
    # 计算指数
    sft <- pickSoftThreshold(exp,
        powerVector = powers,
        networkType = net_type, RsquaredCut = Rsquared_cut, verbose = 5
    )

    # 无标度拓扑拟合指数图，用来选择软阈值，我们一般选择第一个达到Rsquared_cut的power值
    pdf(file = str_glue("{outdir}/Figure_SoftPower.pdf"), width = width, height = height)
    par(mfrow = c(1, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
        main = paste("Scale independence")
    )
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = powers, cex = 0.9, col = "red"
    )
    abline(h = Rsquared_cut, col = "red")
    # Soft threshold与平均连通性图
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
        xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
        main = paste("Mean connectivity")
    )
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
    dev.off()

    # 软阈值筛选
    power <- sft$powerEstimate
    # 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使无标度网络图谱结构R^2达到0.8，
    # 平均连接度较高如在100以上，可能是由于部分样品与其他样品差别太大。
    # 这可能由批次效应、样品异质性或实验条件对表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
    # 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
    if (is.na(power)) {
        message("没有符合条件的软阈值，使用经验power值")
        power <- ifelse(nSamples < 20, ifelse(net_type == "unsigned", 9, 18),
            ifelse(nSamples < 30, ifelse(net_type == "unsigned", 8, 16),
                ifelse(nSamples < 40, ifelse(net_type == "unsigned", 7, 14),
                    ifelse(net_type == "unsigned", 6, 12)
                )
            )
        )
    }
    message(paste0("The pick power: ", power))
    return(list(power = power, sft = sft))
}