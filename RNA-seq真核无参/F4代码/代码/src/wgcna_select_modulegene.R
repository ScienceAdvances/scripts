#' @TODO 根据选定参数筛选WGCNA网络模块基因（核心基因）
#' @title ### 根据选定参数筛选WGCNA网络模块基因（核心基因）
#' @param merged_infor 数据框，wgcna_built_net函数返回的merged_infor，包含基因、模块、以及基因与各模块的MM值、基因与各性状的GS值
#' head(merged_infor)
#     gene module        green       brown     yellow        blue   turquoise Age
#    A1BG   grey -0.359850484 -0.16933946 -0.2263151 -0.45013160 -0.24286994 0.00755854
#    A1CF   grey -0.459555237 -0.17735928 -0.2858948 -0.03973337  0.29127338 -0.06056932
#     A2M   grey -0.052012978  0.23300595  0.1944972 -0.17055515 -0.13229127 -0.12443514
#  A4GALT   grey  0.003539036  0.49099685  0.3131406  0.09260547  0.17013676 0.02833113
#   AADAC   grey -0.374770662 -0.13515580 -0.2620106 -0.34717381 -0.05059942 0.02833113
#   AADAT   grey -0.429576937 -0.09999154 -0.1575306 -0.09508979  0.08668811 0.02833113
#' @param pheno_module_list 列表，包含挑选的性状以及对应的模块，module和pheno的长度要相同
#' pheno_module_list <- list(module = c("green", "brown"), pheno = c("Age", "Age"))
#' @param GS 数值，基因与性状的GS筛选阈值，默认为0不筛选
#' @param MM 数值，基因与模块的MM筛选阈值，默认为0不筛选（GS和MM都不筛选的话，默认为模块全部基因）
#' @param outdir 文件输出路径 
#' @examples wgcna_select_modulegene(merged_infor = wgcna_net$merged_infor, pheno_module_list = list(module = c("green", "brown"), pheno = c("Age", "Age")), outdir = outdir)
#' @author *CY*
#'
wgcna_select_modulegene <- function(merged_infor = NULL, pheno_module_list = NULL ,GS = 0, MM = 0, outdir = getwd()){
    mkdir(outdir)
    res <- lapply(1:length(pheno_module_list$module), function(i) {
        select_module <- pheno_module_list$module[i]
        select_pheno <- pheno_module_list$pheno[i]
        select_data <- merged_infor %>% dplyr::filter(module %in% select_module)
        pdf(file = str_glue("{outdir}/Figure_{select_pheno}_{select_module}_cor.pdf"), width = 6, height = 6)
        verboseScatterplot(abs(select_data[, select_module]),
            abs(select_data[, select_pheno]),
            xlab = paste("Module Membership in", select_module),
            ylab = paste("Gene Significance of ", select_pheno),
            main = paste("Module membership vs Gene Significance\n"),
            cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = select_module, pch = 15
        )
        if (GS > 0) {
            abline(h = GS, col = "red")
        }
        if (MM > 0) {
            abline(v = MM, col = "red")
        }
        dev.off()
        # 筛选模块基因
        select_gene <- select_data %>%
            dplyr::filter(abs(get0(select_module)) >= MM, abs(get0(select_pheno)) >= GS) %>%
            dplyr::pull(gene)
        message(paste0("在GS为",GS,"并且MM为",MM,"时，根据性状",select_pheno,"和模块",select_module,"，共筛选得到",length(select_gene),"个模块基因"))
        return(select_gene)
    })
    names(res) <- pheno_module_list$module
    return(res)
}
