install.packages("immunarch")           # Install the package
library(immunarch); data(immdata)       # Load the package and the test dataset
p=repOverlap(immdata_mixcr$data) %>% vis()      # Compute and visualise the most important statistics:
using::gs(p,outdir="Result/immunarch",name="1",w=12,h=12)
hue=using::hue("NPG")
remove.packages("immunarch")
devtools::install_github("immunomind/immunarch")
BiocManager::install("immunomind/immunarch")
packageVersion("immunarch")

library(immunarch)
using::using(tidyverse)
immdata <- immunarch::repLoad("Result/immunarch/TRUST4")

# public clonotypes, gene usage, sample diversity
geneUsage(immdata$data[[1]]) %>% vis() %>%
    using::gs(outdir="Result/immunarch",name="1",w=12)

# Group samples
repDiversity(immdata$data) %>% vis(.by = "Group", .meta = immdata$meta) %>%
    using::gs(outdir="Result/immunarch",name="2",w=12)

# repExplore
repExplore(immdata$data, .method = "volume") %>% vis() %>%
    using::gs(outdir="Result/immunarch",name="number of unique clonotypes",w=12)

repExplore(immdata$data, .method = "volume") %>% vis(.by = "Group", .meta = immdata$meta) %>%
    using::gs(outdir="Result/immunarch",name="number of unique clonotypes-Group",w=12)

repExplore(immdata$data, .method = "count") %>% vis() %>%
    using::gs(outdir="Result/immunarch",name="distribution of clonotype abundances",w=12)

repExplore(immdata$data, .method = "count") %>% vis(.by = "Group", .meta = immdata$meta) %>%
    using::gs(outdir="Result/immunarch",name="distribution of clonotype abundances-Group",w=12)

repExplore(immdata$data, .method = "len") %>% vis() %>%
    using::gs(outdir="Result/immunarch",name="distribution of CDR3 sequence lengths",w=12)

repExplore(immdata$data, .method = "len") %>% vis(.by = "Group", .meta = immdata$meta) %>%
    using::gs(outdir="Result/immunarch",name="distribution of CDR3 sequence lengths-Group",w=12)


imm_hom <- repClonality(immdata$data,
  .method = "homeo",
  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
fwrite(as.data.table(imm_hom,keep.rownames="Sample"),"Result/immunarch/imm_hom.csv")


repExplore(immdata$data, .method = "len")  %>% fwrite("Result/immunarch/cdr3_len.csv")
repExplore(immdata$data, .method = "volume")  %>% fwrite("Result/immunarch/number of unique clonotypes.csv")
repExplore(immdata$data, .method = "count") %>% fwrite("Result/immunarch/distribution of clonotype abundances.csv")

p=vis(imm_hom, .by ="Group", .meta = immdata$meta) + scale_fill_manual(values=hue)
using::gs(p,outdir="Result/immunarch",name="Clonality-Group",w=12)

# geneUsage HomoSapiens.TRBJ hs.trbv
imm_gu <- geneUsage(immdata$data, "hs.trbv")
p=vis(imm_gu, .by = "Group", .meta = immdata$meta, .plot = "box") + scale_fill_manual(values=hue)
using::gs(p,outdir="Result/immunarch",name="geneUsage-trbv-Group",w=18)

class(imm_gu)