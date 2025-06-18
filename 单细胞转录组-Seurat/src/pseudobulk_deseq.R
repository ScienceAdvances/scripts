pseudobulk_deseq <- function(
    object,
    outdir,
    method = c('average','aggregate')[1],
    sample = "Sample",
    cluster = "CellType",
    condition = "Condition",
    condition1 = "Case",
    condition2 = "Control",
    assay = "RNA") {
    using::mkdir(outdir)
    bls <- object@meta.data %>% dplyr::pull(condition) %in% c(condition1, condition2)
    object %<>% subset(cells = rownames(object@meta.data)[bls])

    pb_object <- Seurat:::PseudobulkExpression(
        object = object,
        assays = assay,
        layer = "counts",
        return.seurat = TRUE,
        group.by = c(condition, sample, cluster),
        method = method[1],
        verbose = TRUE
    )
    saveRDS(pb_object,file=file.path(outdir,'pb_object.rds'))
    strs <- rownames(pb_object@meta.data) %>% stringr::str_split("_", simplify = TRUE)
    
    clusters <- object@meta.data %>% dplyr::pull(cluster)
    cell <- ifelse(ncol(strs)<3,clusters,strs[, 3])
    conditions <- object@meta.data %>% dplyr::pull(condition)
    pb_object$celltype_condition <- paste(cell, strs[, 1], sep = "___")
    Idents(pb_object) <- "celltype_condition"
    pb_object[['RNA']]@layers$counts %<>% apply(c(1,2),round) %>% Matrix::Matrix()
    # x='Astrocytes___Sevoflurane';y='Astrocytes___Control'

    dea_list <- purrr::map2(
        paste(unique(clusters), condition1, sep = "___"), paste(unique(clusters), condition2, sep = "___"),
        function(x, y) {
            tryCatch(
                {
                    tb <- Seurat::FindMarkers(
                        object = pb_object,
                        slot='counts',
                        ident.1 = x,
                        ident.2 = y,
                        logfc.threshold = 0,
                        min.pct = 0,
                        min.cells.feature = 0,
                        min.cells.group = 0,
                        test.use = "DESeq2",
                        pseudocount.use=1
                    )
                    tb %<>% data.table::as.data.table(keep.rownames = "Feature")
                    data.table::setnames(tb, c("Feature", "Pvalue", "LFC", "Pct1", "Pct2", "Padj"))
                    tb %<>% dplyr::arrange(dplyr::desc(LFC)) %>% tidyr::drop_na(LFC,Padj) %>% dplyr::select(-Pct1,-Pct2)
                    return(tb)
                },
                error = function(e) {
                    message('-----------------erorr--------------')
                    return(NULL)
                }
            )
        }
    )
    names(dea_list) <- unique(clusters)
    bls <- map_lgl(dea_list, ~ is.null(.x))
    dea_list <- dea_list[!bls]
    wb <- openxlsx::createWorkbook()
    purrr::walk2(
        dea_list, names(dea_list),
        function(x, y) {
            openxlsx::addWorksheet(wb, sheetName = y)
            openxlsx::writeData(wb, sheet = y, x = x)
        }
    )
    openxlsx::saveWorkbook(wb, base::file.path(outdir, str_glue("{condition1}_vs_{condition2}.xlsx")), overwrite = TRUE)
    return(dea_list)
}
