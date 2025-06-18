library(optparse)

option_list <- list(
  make_option(c("-r", "--rds"), help = "", type = "character", default = "All_sample_combined.rds"),
  make_option(c("-c", "--anno"), help = "cluster_file", type = "character", default = "cluster.txt"),
  make_option(c("-o", "--outdir"), help = "outdir", type = "character", default = "output dir"),
  make_option(c("-t", "--species"), help = "type eg. hsa mmu ..", type = "character", default = "hsa"),
  make_option(c("-p", "--cmpfile"), help = "compare file", type = "character", default = "NULL"),
  make_option(c("-m", "--marker"), help = "each CellType marker for plot", type = "character", default = "NULL"),
  make_option(c("-a", "--avg_log2FC"), help = "threshold for group compare foldchange", type = "double", default = 0.1),
  make_option(c("-s", "--script_dir"), help = "script_dir", type = "character", default = "/home/data/wd/Vik/CODEHUB/SingleCell")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
using::using(using, tidyverse, Seurat, magrittr, patchwork, data.table, cowplot)

# opt <- list(
#   rds = "Result/02_Cluster/Seuart.rds",
#   anno = "Result/04_Annotation/annotation.tsv",
#   avg_log2FC = 0.1,
#   outdir = "Result/04_Annotation",
#   type = "hsa",
#   cmpfile = "NULL",
#   script_dir = "/home/zktbk0/Alex/CODEHUB/SingleCell"
# )

rds <- opt$rds
anno <- opt$anno
outdir <- opt$outdir
species <- opt$species
cmpfile <- opt$cmpfile
marker <- opt$marker
avg_log2FC <- opt$avg_log2FC
script_dir <- opt$script_dir

source('src/config.R')
hue <- hue5

using::mkdir(outdir)

message("============= load anno data =============")
anno <- data.table::fread(opt$anno, sep = "\t", header = T, col.names = c("Cluster", "Anno", "CellTypeDescription", "Marker"))
seurat_anno <- anno %>%
  dplyr::as_tibble() %>%
  tidyr::separate_rows(Cluster, sep = ",")
if (length(seurat_anno$Cluster) != length(unique(seurat_anno$Cluster))) {
  return(message("错误，cluster重复"))
}

new.cluster.ids <- as.vector(seurat_anno$Anno)
names(new.cluster.ids) <- seurat_anno$Cluster


print(seurat_anno$Cluster)

message("============= load seurat data =============")
object <- base::readRDS(opt$rds)
object$Sample %<>% as.factor()
SeuratObject::Idents(object) <- "seurat_clusters"

print(unique(object@active.ident))
object <- subset(object, idents = as.vector(seurat_anno$Cluster))

object <- SeuratObject::RenameIdents(object, new.cluster.ids)
object$CellType <- SeuratObject::Idents(object)

DefaultAssay(object) <- "RNA"
Idents(object) <- object$CellType
saveRDS(object, file = file.path(outdir, "seurat.rds"))


sample_list <- unique(object$Sample)
p1 <- plot.clusters.group(data = object, clusters = "CellType", xlab = "Cluster number", log = TRUE, group = "Group", legend.title = "Group", widths = c(3, 1), hue = hue)
using::gs(p1,outdir = outdir,name = 'Cluster_Group_percent', w = 7, h = 6)

p2 <- plot.clusters.group(data = object, clusters = "CellType", xlab = "Cluster number", log = TRUE, group = "Sample", legend.title = "Sample", widths = c(3, 1), hue = hue)
using::gs(p2,outdir = outdir,name = 'Cluster_Sample_percent', w = 7, h = 6)



# groupDiffAuto(object,opt$outdir,"CellType",opt$type)
# if(opt$cmpfile!= "NULL"){
#     groupDiffSpeci(object,opt$outdir,"CellType",opt$cmpfile,opt$type,opt$avg_log2FC)
# }else{
#     groupDiffAuto(object,opt$outdir,"CellType",opt$type,opt$avg_log2FC)
# }
# cmd <- paste0("cp ",opt$cluster," ",outdir)
# system(cmd)
write.table(table(object$Sample, object$CellType), file = paste(outdir, "Cluster_sample_percent.csv", sep = "/"), sep = ",", quote = FALSE, col.names = NA)
# system(paste0("cp /PERSONALBIO/work/singlecell/s00/software/script/README/添加细胞标签结果反馈说明.p* ",outdir))

# marker图
marker_dotplot(object, marker_anno = anno, outdir = outdir, hue = hue, value_colors = "RdBu",w=12,h=5)

for (i in seq(nrow(anno))) {
  # i=1
  name <- anno$Anno[i]
  gene <- anno$Marker[i] %>%
    str_split(",") %>%
    .[[1]] %>% intersect(rownames(object))
  plots <- map(gene, ~ scCustomize::FeaturePlot_scCustom(
    seurat_object = object,
    colors_use = colorRampPalette(c("#3288BD", "white", "#D53E4F"))(50),
    features = .x
  ) & NoAxes())
  p <- wrap_plots(plots, ncol = 3)
  using::gs(p,outdir = file.path(outdir,"CellMakerUMAP"),name = name, w = 16, h = 7 * ceiling(length(gene) / 3))

}
markers <- scdea(object, outdir = paste0(outdir, "/Diff_Cluster"), hue = hue)
scplot(object, bunch = "CellType", hue = hue, outdir = outdir, marker = markers)
