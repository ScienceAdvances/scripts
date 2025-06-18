library(jsonlite)
library(purrr)
library(RCurl)
library(tidyverse)
library(data.table)
install.packages('ontologyIndex')
# GO功能文件整理
# 读取eggnog注释结果
eggnog <- data.table::fread("orgdb/Z100.emapper.annotations",skip = 4,sep='\t')

# colnames(eggnog)
eggnog$COG_category
# 构建gene id to gene description
gene_info <- eggnog %>%
    dplyr::rename("query_name" = "#query") %>%
    dplyr::select(
        GID = query_name,
        GENENAME = seed_ortholog
    ) %>%
    tidyr::drop_na()

gene2go <- eggnog %>%
    dplyr::rename("query_name" = "#query") %>%
    dplyr::select(GID = query_name, GO = GOs) %>%
    tidyr::separate_rows(GO, sep = ",", convert = F) %>%
    tidyr::drop_na(GO) %>%
    dplyr::filter(str_detect(GO, "^GO")) %>%
    dplyr::mutate(EVIDENCE = "IEA")
head(as.data.frame(gene2go))
gene2go$EVIDENCE %>% table()
gene2go %>% mutate(term=GO,gene=str_remove(GID,'_i[0-9]+.p[0-9]+')) %>%  dplyr::select(term,gene) %>% fwrite('orgdb/go_gene2term.tsv',sep='\t')


go <- ontologyIndex::get_OBO(
  'orgdb/go-basic.obo',
  propagate_relationships = "is_a",
  extract_tags = "everything",
  merge_equivalent_terms = TRUE
)
a <- as.data.frame(go)
fwrite(a,'orgdb/goid_info.tsv',sep='\t')
b <- a %>% select(id,name,namespace)
c <- gene2go %>% mutate(term=GO,gene=str_remove(GID,'_i[0-9]+.p[0-9]+')) %>%
    dplyr::select(term,gene) %>% 
    merge(b,by.x='term',by.y='id')
c$namespace %>% table() # biological_process cellular_component molecular_function
c %>% dplyr::filter(namespace=='biological_process') %>% select(term,gene) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gobp_term2gene.tsv',sep='\t')
c %>% dplyr::filter(namespace=='cellular_component') %>% select(term,gene) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gocc_term2gene.tsv',sep='\t')
c %>% dplyr::filter(namespace=='molecular_function') %>% select(term,gene) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gomf_term2gene.tsv',sep='\t')

c %>% dplyr::filter(namespace=='biological_process') %>% select(term,name) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gobp_term2name.tsv',sep='\t')
c %>% dplyr::filter(namespace=='cellular_component') %>% select(term,name) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gocc_term2name.tsv',sep='\t')
c %>% dplyr::filter(namespace=='molecular_function') %>% select(term,name) %>% distinct(.keep_all = TRUE) %>% fwrite('orgdb/gomf_term2name.tsv',sep='\t')

# KEGG文件整理
# 构建 gene id to kegg pathway关系
# https://www.genome.jp/kegg-bin/get_htext?ko00001
# 整理出KEGG中的所有注释信息，以及大类层级信息

pathway2name <-
    tibble(
        Pathway = character(),
        Name = character(),
        Pathway_subclass = character(),
        Pathway_class = character()
    )

ko2pathway <- tibble(Ko = character(), Pathway = character())

kegg <- fromJSON("orgdb/ko00001.json")
for (a in seq_along(kegg[["children"]][["children"]])) {
    # a=1
    A <- kegg[["children"]][["name"]][[a]]

    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        # b=2
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]

        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
            # c=3
            pathway_info <-
                kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]

            pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
            pathway_name <-
                str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
            Pathway_subclass <- str_remove(B, pattern = "^\\d+ ")
            Pathway_class <- str_remove(A, pattern = "^\\d+ ")
            pathway2name <- rbind(
                pathway2name,
                tibble(
                    Pathway = pathway_id,
                    Name = pathway_name,
                    Pathway_subclass,
                    Pathway_class
                )
            )
            kos_info <-
                kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
            kos <- str_match(kos_info, "K[0-9]*")[, 1]
            ko2pathway <-
                rbind(ko2pathway, tibble(
                    Ko = kos,
                    Pathway = rep(pathway_id, length(kos))
                ))
        }
    }
}

fwrite(pathway2name,'orgdb/pathway2name.tsv',sep='\t')
fwrite(ko2pathway,'orgdb/ko2pathway.tsv',sep='\t')

# 因为是植物相关的研究这里我们将人类相关的通路去除
pathway2name <- fread('orgdb/pathway2name.tsv')%>%
    dplyr::filter(!Pathway_class %in% c("Human Diseases", "Not Included in Pathway or Brite"))


koterms <- eggnog %>% 
    dplyr::select(GID = `#query`, Ko=KEGG_ko)%>% 
    drop_na() %>% 
    dplyr::filter(str_detect(Ko,"ko"))  %>% 
    mutate(Ko=str_remove_all(Ko,"ko:")) %>% 
    tidyr::separate_rows(Ko, sep = ",", convert = F)
colnames(ko2pathway)=c("Ko",'Pathway')
gene2pathway <- koterms %>% left_join(ko2pathway, by = "Ko") %>%dplyr::select(GID, Pathway) %>%drop_na() #合并koterms和ko2pathway到gene2pathway，将基因与pathway的对应关系整理出来
gene2pathway %>% dplyr::filter(Pathway%in%pathway2name$Pathway) %>% mutate(term=Pathway,gene=str_remove(GID,'_i[0-9]+.p[0-9]+')) %>%  dplyr::select(term,gene) %>% fwrite('orgdb/kegg_term2gene.tsv',sep='\t')
pathway2name %>% mutate(term=Pathway,name=str_remove(Name,' \\[BR:ko.+\\]')) %>%  dplyr::select(term,name) %>% fwrite('orgdb/kegg_term2name.tsv',sep='\t')

source("./enrich.R")
for(i in list.files('DESeq2_results','*DE_results$')){
    df=fread(file.path('DESeq2_results',i)) %>% drop_na(padj,log2FoldChange)
    suffix=str_remove_all(i,'(salmon_out.gene.counts.matrix.)|(.DESeq2.DE_results)')
    df %>% dplyr::filter(padj<0.05,log2FoldChange>1) %>% pull(1) %>% enrich(suffix=str_glue("_UP_{suffix}"))
    df %>% dplyr::filter(padj<0.05,log2FoldChange < -1) %>% pull(1) %>% enrich(suffix=str_glue("_DOWN_{suffix}"))
    df %>% dplyr::filter(padj<0.05,abs(log2FoldChange)>1) %>% pull(1) %>% enrich(suffix=str_glue("_ALL_{suffix}"))
}

i='salmon_out.gene.counts.matrix.12h_Light_Chong_vs_0h_chong.DESeq2.DE_results'
df=fread(file.path('DESeq2_results',i)) %>% drop_na(padj,log2FoldChange)
suffix=str_remove_all(i,'(salmon_out.gene.counts.matrix.)|(.DESeq2.DE_results)')
df %>% dplyr::filter(padj<0.05,log2FoldChange>1) %>% pull(1) %>% enrich(suffix=str_glue("_UP_{suffix}"))
df %>% dplyr::filter(padj<0.05,log2FoldChange < -1) %>% pull(1) %>% enrich(suffix=str_glue("_DOWN_{suffix}"))
df %>% dplyr::filter(padj<0.05,abs(log2FoldChange)>1) %>% pull(1) %>% enrich(suffix=str_glue("_ALL_{suffix}"))


