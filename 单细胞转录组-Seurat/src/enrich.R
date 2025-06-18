#' @TODO 富集分析
#' @title ### 功能富集分析
#' @description 调用了子函数`bubble_plot`、`swr`
#' @param suffix 基因类型，用于输出文件的命名
#' @param outdir 结果输出路径
#' @param genelist 基因列表
#' @param hue 字符串向量，代表颜色
#' @return *list*
#' @examples fun_res <- enrich(suffix = "X", genelist = signature, outdir = "./")
enrich <- function(
    genelist,
    suffix = "",
    outdir = ".",
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#6699A1", "#A73D7C")) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(DOSE))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(clusterProfiler))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    suppressPackageStartupMessages(library(xlsx))
    suppressPackageStartupMessages(library(magrittr))

    geneid <- clusterProfiler::bitr(genelist, fromType = fromType, toType = "ENTREZID", OrgDb = orgdb)$ENTREZID
    kegg <- clusterProfiler::enrichKEGG(gene = geneid, organism = species, pvalueCutoff = 1)
    kegg %<>% clusterProfiler::setReadable(OrgDb = orgdb, keyType = "ENTREZID") %>% as.data.frame()
    kegg$Description %<>% stringr::str_remove(pattern = " - .+ \\(.+\\)")
    gobp <- clusterProfiler::enrichGO(gene = geneid, OrgDb = orgdb, ont = "BP", pAdjustMethod = "fdr", readable = TRUE)
    gomf <- clusterProfiler::enrichGO(gene = geneid, OrgDb = orgdb, ont = "MF", pAdjustMethod = "fdr", readable = TRUE)
    gocc <- clusterProfiler::enrichGO(gene = geneid, OrgDb = orgdb, ont = "CC", pAdjustMethod = "fdr", readable = TRUE)

    if (simplify) {
        suppressPackageStartupMessages(library(simplifyEnrichment))
        # 计算相似性矩阵
        mat <- simplifyEnrichment::GO_similarity(gobp@result$ID, ont = "BP", db = "org.Hs.eg.db")
        # 聚类并画图
        pdf(stringr::str_glue("{outdir}/Figure_GO_simplify_{suffix}.pdf"))
        df <- simplifyEnrichment::simplifyGO(mat, method = "binary_cut", plot = TRUE)
        data.table::fwrite(df, file = stringr::str_glue("{outdir}/GO_simplify_{suffix}.tsv"), sep = "\t")
        dev.off()
    }


    ## 画图
    if (nrow(kegg) > 0 && !is.null(kegg)) {
        xlsx::write.xlsx2(kegg, file = stringr::str_glue("{outdir}/Table_KEGG_{suffix}.xlsx"), row.names = FALSE)
        kegg %<>% dplyr::mutate(ID = swr(Description, 50))
        p1 <- bubble_plot(frame = kegg, title = "KEGG Pathway", hue = hue)
        gs(p = p1, outdir = outdir, name = paste0("KEGG_enrich_", suffix), w = 9, h = 6)
    }
    if (nrow(gobp) > 0 && !is.null(gobp)) {
        gobp <- as.data.frame(gobp)
        xlsx::write.xlsx2(gobp, file = stringr::str_glue("{outdir}/Table_GOBP_{suffix}.xlsx"), row.names = FALSE)
        gobp %<>% dplyr::mutate(ID = swr(Description, 50))
        p2 <- bubble_plot(frame = gobp, title = "GO Biological Process", hue = hue)
        gs(p = p2, outdir = outdir, name = paste0("GOBP_enrich_", suffix), w = 9, h = 6)
    }
    if (nrow(gomf) > 0 && !is.null(gomf)) {
        gomf <- as.data.frame(gomf)
        xlsx::write.xlsx2(gomf, file = stringr::str_glue("{outdir}/Table_GOMF_{suffix}.xlsx"), row.names = FALSE)
        gomf %<>% dplyr::mutate(ID = swr(Description, 50))
        p3 <- bubble_plot(frame = gomf, title = "GO Molecular Function", hue = hue)
        gs(p = p3, outdir = outdir, name = paste0("_GOMF_enrich_", suffix), w = 9, h = 6)
    }
    if (nrow(gocc) > 0 && !is.null(gocc)) {
        gocc <- as.data.frame(gocc)
        xlsx::write.xlsx2(gocc, file = stringr::str_glue("{outdir}/Table_GOCC_{suffix}.xlsx"), row.names = FALSE)
        gocc %<>% dplyr::mutate(ID = swr(Description, 50))
        p4 <- bubble_plot(frame = gocc, title = "GO Cellular Component", hue = hue)
        gs(p = p4, outdir = outdir, name = paste0("GOCC_enrich_", suffix), w = 9, h = 6)
    }
    ## 合图
    if (exists("p1") & exists("p2") & exists("p3") & exists("p4")) {
        p <- cowplot::plot_grid(p1, p2, p3, p4, labels = NA, ncol = 2, label_size = 24, scale = c(1, 1, 1, 1))
        gs(p = p, outdir = outdir, name = paste0("Enrich_", suffix), w = 18, h = 12)
    }
    return(list(KEGG = kegg, GOBP = gobp, GOCC = gocc, GOMF = gomf))
}

bubble_plot <- function(frame, title, hue) {
    frame %<>% dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
        dplyr::mutate(Padjust = as.numeric(p.adjust))
    frame %<>% dplyr::mutate(FoldEnrichment = DOSE::parse_ratio(GeneRatio) / DOSE::parse_ratio(BgRatio))
    p <- ggplot2::ggplot(data = frame, aes(FoldEnrichment, y = reorder(ID, -log10(Padjust)))) +
        ggplot2::geom_point(aes(colour = -log10(Padjust), size = Count)) +
        scale_size_continuous(range = c(2, 6)) +
        scale_colour_gradient(low = hue[1], high = hue[2]) +
        ggtitle(title) +
        ylab("") +
        xlab("Rich factor") +
        theme_bw() +
        theme(
            plot.title = element_text(vjust = 1, hjust = 0.5),
            legend.key = element_blank(),
            title = element_text(face = "bold", size = 15),
            axis.text.x = element_text(face = "bold", size = 10, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            legend.title = element_text(face = "bold", size = 8),
            legend.text = element_text(size = 8)
        )
    return(p)
}
