using::using(using,ChIPseeker,GenomicFeatures,tidyverse,ggupset)

name="KO"
outdir="Result/05.Annotation"
tssRegion=c(-3000, 3000)
colors=NA
peak_format="narrow"
meta <- Sys.getenv("META")

cat = data.table::fread(meta, sep='\t') %>% dplyr::filter(Group!="input") %>% dplyr::pull(Group) %>% unique()
message("groups: ", cat)

genome <- GenomicFeatures::makeTxDbFromGFF('DB/mmu.v37.gencode.pri.gtf',format="gtf")


for(name in cat){
    for(peak_format in c("narrow","broad")){
    peak <- ChIPseeker::readPeakFile(glue::glue('Result/03.Peak_Calling/{name}_peaks_{peak_format}Peak_IDR.bed'))
    peakAnno <- ChIPseeker::annotatePeak(peak, tssRegion=tssRegion, TxDb=genome)
    as.data.frame(peakAnno) %>% data.table::fwrite(glue::glue("Result/05.Annotation/{name}_peakAnnotation_{peak_format}.xls"),sep='\t')

    p_covplot <- ChIPseeker::covplot(peak, weightCol="V5")
    using::gs(p_covplot,name=glue::glue("{name}_covplot_{peak_format}"),outdir = outdir)

    p_upsetplot <- ChIPseeker::upsetplot(peakAnno, vennpie=TRUE)
    using::gs(p_upsetplot,name=glue::glue("{name}_upsetplot_{peak_format}"),outdir = outdir)

    p_plotDistToTSS <- ChIPseeker::plotDistToTSS(peakAnno,title="Distribution of TF-binding loci relative to TSS")
    using::gs(p_plotDistToTSS,name=glue::glue("{name}_plotDistToTSS_{peak_format}"),outdir = outdir)

    p_AnnoBar <- ChIPseeker::plotAnnoBar(peakAnno)
    using::gs(p_AnnoBar,name=glue::glue("{name}_AnnoBar_{peak_format}"),outdir = outdir)


    pdf(glue::glue("{outdir}/{name}_AnnoPie_{peak_format}.pdf"),width = 9,height=6)
    ChIPseeker::plotAnnoPie(peakAnno,col=NA)
    dev.off()

    pdf(glue::glue("{outdir}/{name}_VennPie_{peak_format}.pdf"),width = 9,height=6)
    ChIPseeker::vennpie(peakAnno,col=NA)
    dev.off()
}}
