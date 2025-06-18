Rscript src/run_subcluster.R \
--rds Result/04_Annotation/seurat.rds \
--resolution 0.6 \
--outdir Result/05_Subset/MPs/Cluster \
--ident CellType \
--subset MPs

Rscript src/addlable.R \
--rds Result/05_Subset/MPs/Cluster/Seurat.rds \
--anno Result/05_Subset/MPs/Annotation/annotation.tsv \
--outdir Result/05_Subset/MPs/Annotation \
--species hsa
