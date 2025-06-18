message(crayon::bgCyan("=============== Package Version ==============="))
message("Seurat: ", packageVersion("Seurat"))
message("SeuratObject: ", packageVersion("SeuratObject"))
using::using(ComplexHeatmap,glue,GseaVis)

# =============== 导入数据 ===============
WORKDIR <- "/home/alex/Seurat"
source("src/config.R")
outdir="Result/06.Pseudobulk"
srt <- readRDS('Result/06_Reflection/Seurat.rds')

pseudobulk_deseq(
    object = srt,
    outdir = outdir,
    method = c('average','aggregate')[2],
    sample = "Sample",
    cluster = "CellType",
    condition = "Group",
    condition1 = "Case",
    condition2 = "Control",
    assay = "RNA"
)



deg_table <- readxl::read_xlsx("S147-2/NUPR1_Pos_vs_NUPR1_Neg.xlsx")
deg_table %<>% dplyr::filter(Pvalue < 0.05, abs(LFC) > 0.585) %>% arrange(desc(abs(LFC)),Pvalue)
gene <- deg_table$Feature

# 热图
pb_object <- Seurat:::PseudobulkExpression(
    object = srt,
    assays = "RNA",
    layer = "counts",
    return.seurat = TRUE,
    group.by = c("Sample", "Group", "CellType"),
    method = "aggregate",
    verbose = TRUE
)
suppressPackageStartupMessages(library())
plotdata <- SeuratObject::GetAssayData(pb_object, layer = "scale.data", assay = "RNA")
colnames(plotdata) %<>% str_remove_all("(^g)|(Normal)|(Tumor)|(_Macrophages)")

pdata <- data.frame(Group = str_split(colnames(pb_object), "_", simplify = T)[, 1], row.names = colnames(plotdata))
pdata %<>% arrange(Group)
plotdata <- plotdata[gene, rownames(pdata)] 

deg_table <- readxl::read_xlsx(file.path(outdir,"Case_vs_Control.xlsx"))
deg_table %<>% dplyr::filter(Pvalue < 0.05) %>% arrange(desc(abs(LFC)))

g1 <- deg_table %>% slice_min(order_by=LFC,n = 5) %>% pull(1)
g2 <- deg_table %>% slice_max(order_by=LFC,n = 5) %>%    pull(1)
gene <- c(g1,g2)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = pdata$Group,
    col = list(Group = c(
        "Case" = hue[3],
        "Control" = hue[1]
    )),
    gp = grid::gpar(col = NA)
)

p <- ComplexHeatmap::Heatmap(
    as.matrix(plotdata)[gene,],
    top_annotation = top_annotation,
    column_order = rownames(pdata),
    col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
    name = "zscore",
    show_row_names = T,
    show_column_names = T,
    column_title_gp = grid::gpar(fill = c(case_label = hue5[3], control_label = hue5[1]))
)
using::gs(p = p, name = "heatmap", h = 6, w = 12, outdir = outdir)

# enrich
dplyr::filter(deg_table, abs(LFC) > 0.585, Padj < 0.05) %>%
    pull(1) %>%
    enrich(json = "src/Enrich/enrich.json", species = "hsa", suffix = glue("All"), outdir = outdir)

# gsea
kegg <- read.gmt("src/Enrich/KEGG_2021_Human.gmt")
gobp <- read.gmt("src/Enrich/GO_Biological_Process_2023.gmt")
gocc <- read.gmt("src/Enrich/GO_Cellular_Component_2023.gmt")
gomf <- read.gmt("src/Enrich/GO_Molecular_Function_2023.gmt")
wikipathway <- read.gmt("src/Enrich/WikiPathways_2024_Human.gmt")
reactome <- read.gmt("src/Enrich/Reactome_Pathways_2024.gmt")
geneset <- list(kegg, gobp, gocc, gomf, wikipathway, reactome)
names(geneset) <- c("kegg", "gobp", "gocc", "gomf", "wikipathway", "reactome")

main <- function(iname, geneset, prefix) {
    # sheets <- readxl::excel_sheets(iname)
    plotx(deg_table = fread(iname), geneset = geneset, prefix = prefix)
}
plotx <- function(deg_table, geneset, prefix) {
    # deg_table <- readxl::read_xlsx(deg_table,sheet = "PP")
    deg_table %<>% as.data.frame() %>% drop_na(LFC) %>% dplyr::distinct(Feature, .keep_all = TRUE)
    genelist <- deg_table$LFC
    names(genelist) <- deg_table$Feature %>% str_to_upper()
    genelist %<>% sort(decreasing = T)
    res <- map2(names(geneset), geneset, function(x, y) {
        # pathways <- map(split(y, y$term), ~ pull(.x, 2))
        a <- clusterProfiler::GSEA(geneList=genelist,TERM2GENE=y,pvalueCutoff=1)
        # a <- fgsea(
        #     pathways = pathways,
        #     stats = genelist,
        #     minSize = 10,
        #     maxSize = 500,
        #     nperm = 1000
        # )
        saveRDS(a,file=base::file.path(outdir, str_glue("{x}.rds")))
        return(as.data.frame(a))
    })
    wb <- openxlsx::createWorkbook()
    purrr::walk2(
        res, names(geneset),
        function(x, y) {
            openxlsx::addWorksheet(wb, sheetName = y)
            openxlsx::writeData(wb, sheet = y, x = x)
        }
    )
    openxlsx::saveWorkbook(wb, base::file.path(outdir, str_glue("{prefix}.xlsx")), overwrite = TRUE)
}
# for (x in list.files("../差异分析结果")) {
#     main(iname = glue::glue("../差异分析结果/{x}"), geneset = geneset, prefix = str_remove(x, ".csv"))
# }
deg_table <- readxl::read_xlsx(file.path(outdir,"Case_vs_Control.xlsx"))
plotx(deg_table = deg_table, geneset = geneset, prefix = "Macrophages")


gobp <- readRDS(base::file.path(outdir, "gobp.rds"))
p <- GseaVis::gseaNb(object=gobp,addPval=T,geneSetID='Response To Reactive Oxygen Species (GO:0000302)')
using::gs(p,outdir = outdir,name='ROS')
