source("/home/data/wd/Vik/CODEHUB/RCode/config.R")
g=fread('F4/WGCNA/HubGene/greenyellow.csv') %>% pull(1)
g=fread('F4/WGCNA/HubGene/salmon.csv') %>% pull(1)
g=fread('F4/WGCNA/HubGene/purple.csv') %>% pull(1)
g=fread('F4/WGCNA/HubGene/yellow.csv') %>% pull(1)
tf=fread('F2/数据/tf_classification.txt',header=FALSE)
tf=tf[V3=='TF',]
tf$Gene <- map_chr(str_split(tf$V1,'_'),~paste0(.x[1:4],collapse = '_'))
var=c('V2','Gene')

for(i in c('brown','green','red','yellow','turquoise')){
    g=fread(glue::glue('F4/WGCNA/HubGene/{i}.csv')) %>% pull(1)
    net <- tf[Gene %chin% g,..var] %>% distinct(.keep_all=TRUE)
    attr <- data.table(node=c(net$V2,net$Gene),attr=c(rep('TF',length(net$V2)), rep('Gene',length(net$Gene)))) %>% 
        distinct(.keep_all=TRUE)
    fwrite(net, glue::glue('F4/Net/{i}N.csv'))
    fwrite(attr,  glue::glue('F4/Net/{i}A.csv'))
}