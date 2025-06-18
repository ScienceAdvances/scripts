mamba env create -f /home/data/wd/Vik/S134/citeseq.yaml
mamba activate citeseq
CITE-seq-Count \
--read1 MGS_R1.fq.gz \
--read2 MGS_R2.fq.gz \
--tags meta/tag.csv \
--cell_barcode_first_base 1 \
--cell_barcode_last_base 16 \
--umi_first_base 17 \
--umi_last_base 28 \
--expected_cells 10000 \
--whitelist meta/3M-february-2018.txt \
--output MGS \
--threads 32

cellranger-8.0.1/cellranger count \
--id PBS \
--libraries meta/PBS.csv \
--feature-ref meta/feature-ref.csv \
--transcriptome ~/APP/refdata-gex-GRCm39-2024-A \
--localmem 100 \
--localcores 36 \
--create-bam true