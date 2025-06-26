mamba env create -f src/chip.yaml
mamba activate chip
Rscript -e 'BiocManager::install("ScienceAdvances/using",update=F,ask=F)'

bowtie2-build --threads 32 mmu.v37.gencode.pri.fna mmu