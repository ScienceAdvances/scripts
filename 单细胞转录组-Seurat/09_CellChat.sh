Rscript src/run_CellChat.R \
--seurat Result/06_Reflection/Seurat.rds \
--filename Case \
--outdir Result/09_CellChat \
--species hsa \
--workdir Result/09_CellChat \
--n_jobs 16 \
--cluster Group \
--celltype CellType \
--conda /home/data/wd/miniforge3/envs/sc \
--used_celltype Result/09_CellChat/used_celltype.txt