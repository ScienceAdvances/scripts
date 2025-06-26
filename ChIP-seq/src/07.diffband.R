library(optparse)
option_list <- list(
    make_option("--specie", type = "character", action = "store", default = "mmu", help = "mmu/hsa"),
    make_option("--peak_format", type = "character", action = "store", default = "narrow", help = "narrow/broad")
)
opt <- parse_args(OptionParser(
    option_list = option_list, add_help_option = TRUE,
    usage = "Usage: %prog [options] \nDescription: CellChat pre!"
))
specie=opt$specie
peak_format=opt$peak_format
message("===================================== specie:", specie)
message("===================================== peak_format:", peak_format)
using::using(tidyverse,ChIPQC,DiffBind,BiocParallel,GenomicRanges,GenomicFeatures )
BiocParallel::register(BiocParallel::MulticoreParam(workers = 32))

meta <- read.table('experiment.tsv',header = T)
meta$bamReads = paste0("Result/02.Alignment/",meta$SampleID,".bam")
meta$bamControl = paste0("Result/02.Alignment/",meta$ControlID,".bam")
meta$Factor = meta$Group

if(peak_format=="narrow"){
    meta$PeakCaller = "narrow"
    meta$Peaks =  glue::glue("Result/03.Peak_Calling/{peak_format}Peak/{meta$SampleID}_peaks.narrowPeak")
}else if(peak_format=="broad"){
    meta$PeakCaller = "macs"
    meta$Peaks =  glue::glue("Result/03.Peak_Calling/{peak_format}Peak/{meta$SampleID}_peaks.xls")
}
print(meta)
txdb=switch(specie,
    mmu=TxDb.Mmusculus.UCSC.mm39.knownGene::TxDb.Mmusculus.UCSC.mm39.knownGene,
    hsa=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
)
All5utrs <- reduce(unique(unlist(fiveUTRsByTranscript(txdb))))
All3utrs <- reduce(unique(unlist(threeUTRsByTranscript(txdb))))
Allcds <- reduce(unique(unlist(cdsBy(txdb,"tx"))))
Allintrons <- reduce(unique(unlist(intronsByTranscript(txdb))))
Alltranscripts <- reduce(unique(transcripts(txdb)))

posAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "+"]
posAllTranscripts <- posAllTranscripts[!(start(posAllTranscripts)-20000 < 0)]
negAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "-"]
chrLimits <- seqlengths(negAllTranscripts)[as.character(seqnames(negAllTranscripts))]      
negAllTranscripts <- negAllTranscripts[!(end(negAllTranscripts)+20000 > chrLimits)]      
Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
Promoters500 <-  reduce(flank(Alltranscripts,500))    
Promoters2000to500 <-  reduce(flank(Promoters500,1500))
LongPromoter20000to2000  <- reduce(flank(Promoters2000to500,18000))

annotation=list(
    version="annotation",
    LongPromoter20000to2000=LongPromoter20000to2000,
    Promoters2000to500=Promoters2000to500,
    Promoters500=Promoters500,
    All5utrs=All5utrs,
    Alltranscripts=Alltranscripts,
    Allcds=Allcds,
    Allintrons=Allintrons,
    All3utrs=All3utrs
)

chipObj = ChIPQC::ChIPQC(meta,annotation=annotation)

dbObj <- DiffBind::dba(chipObj)
dbObj <- DiffBind::dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
outdir="Result/07.DIff"
pdf(file.path(outdir,glue::glue("PCA_{peak_format}.pdf")))
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

pdf(file.path(outdir,glue::glue("correlation_heatmap_{peak_format}.pdf")))
plot(dbObj)
dev.off()

dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)

pdf(file.path(outdir,glue::glue("PCA_{peak_format}_sig.pdf")))
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

pdf(file.path(outdir,glue::glue("venn_{peak_format}.pdf")))
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()

pdf(file.path(outdir,glue::glue("MA_DESEQ2_{peak_format}.pdf")))
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()

pdf(file.path(outdir,glue::glue("boxplot_{peak_format}.pdf")))
pvals <- dba.plotBox(dbObj)
dev.off()

dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)  %>% as.data.frame() %>% 
    writexl::write_xlsx(file.path(outdir,glue::glue("edgeR_{peak_format}.xlsx")))

options(browser = "firefox")
a=ChIPQCreport(chipObj, reportName="ChIP QC report", facetBy="Factor",reportFolder="Result/10.stat/ChIPQCreport")
