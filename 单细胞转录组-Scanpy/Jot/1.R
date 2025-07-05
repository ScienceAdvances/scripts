load('/Users/victor/Downloads/scCustomize-master/data/ensembl_mito_id.rda')
load('/Users/victor/Downloads/scCustomize-master/data/ensembl_ribo_id.rda')
str(ensembl_mito_id,1)
str(ensembl_ribo_id,1)
using(scCustomize,clusterProfiler,org.Hs.eg.db,biomaRt)
methods('Add_Mito_Ribo.Seurat')
?Add_Mito_Ribo.Seurat
getAnywhere(Add_Mito_Ribo.Seurat)
scCustomize:::Add_Mito_Ribo.Seurat
clusterProfiler::bitr(geneID=ensembl_mito_id$Homo_sapiens_mito_ensembl, fromType='ENSEMBL', toType='SYMBOL', OrgDb='org.Hs.eg.db', drop = F)

ensembl = biomaRt::useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl", mirror = "aisa")
new=getBM(attributes=c("ensembl_transcript_id_version","refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), values = ensembl_mito_id$Homo_sapiens_mito_ensembl, mart= ensembl)