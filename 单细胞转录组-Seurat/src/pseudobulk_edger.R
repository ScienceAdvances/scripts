pseudobulk_edger <- function(
    object,
    outdir,
    prefix,
    method = c('average','aggregate')[2],
    sample = "Sample",
    cluster = "CellType",
    used_cell,
    condition = "Condition",
    condition1 = "Case",
    condition2 = "Control",
    lfc_threshold = 0.585,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c("#84bbfb","#fc85b8",  "grey60"), #c("#6699A1", "#A73D7C", "#A6A6A6"),
    hue2 = c("#84bbfb", "white","#fc85b8"),
    assay = "RNA") {
    using::mkdir(outdir)
    cells <- object@meta.data %>% dplyr::filter(!!sym(condition) %in% c(condition1, condition2), !!sym(cluster) %in% used_cell) %>% rownames()
    object %<>% subset(cells = cells)
    
    pb_object <- Seurat:::PseudobulkExpression(
        object = object,
        assays = assay,
        layer = "counts",
        return.seurat = TRUE,
        group.by = c(condition, sample),
        method = method,
        verbose = TRUE
    )

    fdata=SeuratObject::GetAssayData(pb_object,layer='counts') %>% as.data.frame()
    fdata %>% as.data.table(keep.rownames='Feature') %>% fwrite(file.path(outdir,glue::glue('{prefix}_fdata.csv')))
    pdata=object@meta.data %>% 
        rownames_to_column('CELLS') %>% 
        dplyr::mutate(ROWNAMES=paste(!!sym(sample),!!sym(condition),sep='_')) %>% 
        dplyr::select(ROWNAMES,!!sym(condition)) %>% 
        dplyr::distinct() %>% 
        column_to_rownames(var='ROWNAMES')  
    deg_table <- edger(
        outdir=outdir,
        prefix=prefix,
        fdata=fdata,
        pdata=pdata,
        control_label=condition2,
        case_label=condition1,
        lfc_threshold = lfc_threshold,
        pvalue_threshold = pvalue_threshold,
        pvalue = pvalue,
        text_gene = NULL,
        hue1 = hue1,
        hue2 = hue2
    )
}
