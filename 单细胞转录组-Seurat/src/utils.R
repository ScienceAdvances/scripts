create_sc <- function(x, y) {
    colnames(x) <- paste(colnames(x), y, sep = "_")
    sc <- SeuratObject::CreateSeuratObject(
        counts = x,
        min.cells = min_cells,
        min.features = min_features,
        project = y
    )
    return(sc)
}