source("src/config.R")
using(Seurat,optparse, clusterProfiler,reshape2,cowplot,scales,patchwork,paletteer,future, RColorBrewer,  infercnv, tidyverse, future, phylogram, gridExtra, grid, dendextend, ggthemes, miscTools)

option_list <- list(
  make_option(c("-m","--mainrds"), help="大类rds/主rds"),
  make_option(c("-s","--c"),help="subrds,use , split"),
  make_option(c("-o","--outdir"),help="outdir",default="./"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa"),
  make_option(c("-p","--step"),help="need rds or all(rds,umap,marker)",default="all"),
  make_option(c("-x","--script_dir"),help="need rds or all(rds,umap,marker)",default="src")
  )
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mainrds <- opt$mainrds
subrds <- opt$subrds
outdir <- opt$outdir
step <- opt$step
type <- opt$type
script_dir <- opt$script_dir

# mainrds <- "Result/04_Annotation/seurat.rds"
# subrds <- "c"
# outdir <- "Result/reflect"
# step = "all"
# type = "hsa"
# hue = hue3

using::mkdir(outdir)
main <- readRDS(mainrds)
subrds <- strsplit(subrds, split = ",") %>% unlist()

colnames(main) %>% head()

# strs <- str_split(colnames(main),'_',simplify=T)

# main %<>% SeuratObject::RenameCells(new.names=paste(strs[,1],strs[,ncol(strs)],sep='_'))
for (i in subrds){
  scsub <- readRDS(i)
  print(i)
  message(colnames(scsub) %>% head())

  # scsub %<>% SeuratObject::RenameCells(new.names=paste(strs[,1],strs[,ncol(strs)],sep='_'))
  main$Barcode <- colnames(main)
  scsub$Barcode <- colnames(scsub)
  main$CellType %<>% as.character()
  scsub$CellType %<>% as.character()
  meta <- scsub@meta.data
  for (j in unique(scsub$CellType)){
    a <- subset(meta, CellType == j)
    main@meta.data[rownames(a),'CellType'] <- j
    print(unique(main$CellType))
}
}

SeuratObject::Idents(main) <- "CellType"
celltypes <- SeuratObject::Idents(main) %>% levels()
# main %<>% subset(idents = base::setdiff(celltypes,c("Macrophages")))
SeuratObject::Idents(main) <- "CellType"
saveRDS(main, file = base::file.path(outdir,"Seurat.rds"))

# =============== find marker ===============
SeuratObject::DefaultAssay(main) <- "RNA"

using::mkdir(file.path(outdir, "Diff_Cluster"))
markers <- scdea(main, outdir = file.path(outdir, "Diff_Cluster"), hue = hue)
scplot(main, bunch = "CellType", hue = hue, outdir = file.path(outdir, "Diff_Cluster"), marker = markers)

p1 = plot.clusters.group(data = main,group='Group',clusters =  "CellType", 
xlab = "Cluster number", log =TRUE, legend.title = "Group",widths = c(3,1),hue = hue)
using::gs(p1,outdir=outdir,name="Cluster_group_percent",w = 12,h = 12,format=c('png','pdf'))

p2 = plot.clusters.group(data = main,clusters =  "CellType", xlab = "Cluster number", 
log =TRUE, group = "Sample",legend.title = "Sample",widths = c(3,1),hue = hue)
using::gs(p2,outdir=outdir,name="Cluster_sample_percent",w = 20,h = 20)

