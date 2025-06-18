#' @TODO 根据选定参数构建WGCNA网络
#' @title ### 根据选定参数构建WGCNA网络
#' @param exp 表达谱文件,样本在行，基因在列
#' @param outdir 文件输出路径
#' @param power 软阈值，wgcna_picksoftpower函数返回值
#' @param pheno 临床信息表，包含sample列，其它列为分析用的性状列
#' @param maxBlockSize 计算机能处理的最大模块的基因数量 ，默认24000；4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可处理3万个
#' @param cor_type 相关性算法，可选pearson和bicor，默认pearson；后者更能考虑离群点的影响
#' @param TOMType 默认为"unsigned"，计算TOM矩阵时，是否考虑正负相关性；可选"signed"，但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义
#' @param minModuleSize 模块中最少的基因数
#' @param mergeCutHeight 合并模块时的聚类树切割距离，设置合并相似性模块的距离，值越小，就越不容易被合并，保留下来的模块就越多
#' @param saveTOMs 逻辑值，默认TRUE，是否保存TOM矩阵；如果设为TRUE，需要设置saveTOMFileBase参数，提供保存文件名
#' @param deepSplit 模块划分的敏感度，0-4之间的整数，参数调整划分模块的敏感度；值越大越敏感，得到的模块就越多，默认是2
#' @param saveTOMFileBase 保存TOM矩阵的文件名
#' @param loadTOM 逻辑值，是否加载先前跑的TOM矩阵，可节省大量运算时间
#' @return *list*
#' @examples wgcna_built_net(exp = qc_res$use_exp, pheno = qc_res$use_pheno, power = softpower_res$power, outdir = outdir, width = 4)
#' @author *CY*
#'
wgcna_built_net <- function(exp = NULL, pheno = NULL, power = NULL, maxBlockSize = 24000, TOMType = "signed",networkType = "unsigned",
                            cor_type = "pearson", minModuleSize = 20, mergeCutHeight = 0.15,
                            saveTOMs = TRUE, deepSplit = 2, outdir = getwd(), saveTOMFileBase = "TOM", loadTOM = FALSE,
                            pamRespectsDendro = FALSE, reassignThreshold = 0, nThreads = 8,hue2=c("#6699A1",'white', "#A73D7C"),
                            width = 4, height = 12, verbose = 3) {
    sam_order <- intersect(rownames(exp), rownames(pheno))
    exp <- exp[sam_order, ]
    pheno <- pheno[sam_order, , drop = FALSE]
    using::using(WGCNA)
    options(bitmapType = "cairo")
    options(stringsAsFactors = FALSE)
    using::mkdir(outdir)
    # 开启多线程
    message("============ built net ============")
    net <- blockwiseModules(
        datExpr = exp,
        power = power,
        nThreads = nThreads,
        maxBlockSize = maxBlockSize,
        networkType=networkType,
        TOMType = TOMType,
        corType = cor_type,
        randomSeed = 1314520,
        minModuleSize = minModuleSize,
        mergeCutHeight = mergeCutHeight,
        numericLabels = FALSE,
        saveTOMs = saveTOMs,
        deepSplit = deepSplit,
        saveTOMFileBase = file.path(outdir, saveTOMFileBase),
        loadTOM = loadTOM,
        verbose = verbose
    )

    moduleColors <- net$colors

    # 统计模块基因数目
    module_gene_count <- table(moduleColors) %>%
        as.data.frame() %>%
        rename(module = moduleColors, number = Freq)
    write.table(module_gene_count, file = file.path(outdir, "SupplementaryTable_module_gene_count.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    # 各模块基因列表
    module_genes <- moduleColors %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        rename(module = 2)
    write.table(module_genes, file = file.path(outdir, "SupplementaryTable_module_genes.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

    # 层级聚类树展示各个模块
    pdf(file = str_glue("{outdir}/Figure_ClusterDendrogram.pdf"), width = 8, height = 6)
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        groupLabels = c(
            "Module colors",
            "GS.weight"
        ),
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05
    )
    dev.off()

    # 关联表型性状
    # 模块性状相关性图
    robustY <- ifelse(cor_type == "pearson", T, F)
    if (cor_type == "pearson") {
        modTraitCor <- stats::cor(net$MEs, pheno, use = "p", method = "pearson")
        modTraitP <- corPvalueStudent(modTraitCor, nrow(pheno))
    } else {
        modTraitCorP <- bicorAndPvalue(net$MEs, pheno, robustY = robustY)
        modTraitCor <- modTraitCorP$bicor
        modTraitP <- modTraitCorP$p
    }
    write.table(modTraitCor, file = file.path(outdir, "SupplementaryTable_modTraitCor.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    write.table(modTraitP, file = file.path(outdir, "SupplementaryTable_modTraitP.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    # 热图可视化模块与表型数据的相关系数和P值
    textMatrix <- paste(sprintf("%.4f", modTraitCor), "\n(", sprintf("%.4f", modTraitP), ")", sep = "")
    dim(textMatrix) <- dim(modTraitCor)
    pdf(str_glue("{outdir}/ModuleTrait.pdf"), width = width, height =height)
    par(mar = c(5, 6, 3, 3)) # 下左上右
    label <- substring(rownames(modTraitP), 3)
    WGCNA::labeledHeatmap(
        Matrix = modTraitCor,
        xLabels = colnames(pheno),
        yLabels = label,
        cex.lab = 1,
        colors= blueWhiteRed(50),
        # colors=grDevices::colorRampPalette(c('blue','white','red'))(100),
        ySymbols = label,
        colorLabels = FALSE,
        textMatrix = textMatrix, 
        setStdMargins = FALSE,
        cex.text = 0.6,
        zlim = c(-1, 1),
        main = paste("Module Trait Relationships")
    )
    dev.off()

    # 基因与模块的MM值
    geneModuleMembership <- cor(exp, net$MEs, use = "p", method = "pearson")
    colnames(geneModuleMembership) <- substring(colnames(geneModuleMembership), 3)
    write.table(geneModuleMembership, file = file.path(outdir, "SupplementaryTable_geneModuleMembership.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    # 基因与性状的GS值
    geneTraitSignificance <- cor(exp, pheno, use = "p", method = "pearson")
    write.table(geneTraitSignificance, file = file.path(outdir, "SupplementaryTable_geneTraitSignificance.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

    # 针对每个性状,每个模块基因的显著性值
    mkdir(str_glue("{outdir}/ModuleSignificance"))
    for (i in colnames(pheno)) {
        geneSignificance <- abs(geneTraitSignificance[, i, drop = FALSE])
        ModuleSignificance_mean <- tapply(geneSignificance, moduleColors, mean, na.rm = T)
        ModuleSignificance_se <- tapply(geneSignificance, moduleColors, stdErr)
        y_max <- max(ModuleSignificance_mean + 1.96 * ModuleSignificance_se) + 0.05
        pdf(str_glue("{outdir}/ModuleSignificance/Figure_{i}_ModuleSignificance.pdf"), width = 10, height = 6)
        barplot(ModuleSignificance_mean,
            names.arg = names(table(moduleColors)), ylim = c(0, y_max),
            col = names(table(moduleColors)), ylab = "Gene Significance", main = "Gene significance across modules"
        )
        addErrorBars(as.vector(ModuleSignificance_mean), as.vector(1.96 * ModuleSignificance_se), two.side = TRUE)
        dev.off()
        ModuleSignificance <- as.data.frame(ModuleSignificance_mean) %>% arrange(desc(ModuleSignificance_mean))
        write.table(ModuleSignificance, file = str_glue("{outdir}/ModuleSignificance/SupplementaryTable_{i}_ModuleSignificance.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    }

    merged_infor <- module_genes %>%
        merge(., geneModuleMembership %>% as.data.frame() %>% rownames_to_column("gene")) %>%
        merge(., geneTraitSignificance %>% as.data.frame() %>% rownames_to_column("gene"))
    # 返回结果
    return(list(
        net = net, module_genes = module_genes, geneModuleMembership = geneModuleMembership,
        geneTraitSignificance = geneTraitSignificance, merged_infor = merged_infor
    ))
}
