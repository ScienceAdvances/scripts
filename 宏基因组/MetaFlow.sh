#/usr/bin/env bash
export WORKDIR=/home/alex/MetaFlow
cd $WORKDIR
source $WORKDIR/scripts/function.sh
export -f $(declare -F | awk '{print $NF}')
export MetaFlowDB=/mnt/e/Alex/X207/MetaFlowDB
META=/mnt/e/Alex/X207/meta.tsv
export THREADS=45
export bwa_index=/mnt/e/Alex/DATAHUB/index/hsa_bwa/hsa
export bowtie2_index=/mnt/e/Alex/DATAHUB/index/hsa_bowtie2/hsa
megan=/home/data/wd/Vik/APP/megan/tools
mamba activate MetaFlow

mkdir -p tmp Result/01.clean_data/ Result/02.assembly tmp Result/03.gene_prediction Result/04.annotation Result/05.taxonomy Result/06.statistics/fastp_report

# Get the list of samples from the meta.tsv file
mapfile -t samples < <(awk 'NR>1 {print $1}' $META)
export samples
mamba install -c bioconda  -c conda-forge -y fastp
# Run fastp on each sample
print_sample | parallel -j 16 run_fastp

# remove host reads from each sample
print_sample | parallel -j 4 run_bwa
# print_sample | parallel -j 4 run_bowtie2

# co
megahit \
-1 $(printf "Result/01.clean_data/%s_1.fq.gz," "${samples[@]}") \
-2 $(printf "Result/01.clean_data/%s_2.fq.gz," "${samples[@]}") \
--presets meta-large \
--out-dir tmp \
--out-prefix megahit \
--num-cpu-threads 32 \
--memory 0.7 \
--verbose
mv tmp/R/contigs.fa tmp/contigs.fa

# merged mode
print_sample | parallel -j 2 run_megahit
print_sample | parallel -j 8 rename_megahit
cat tmp/*/megahit.contigs.fa > tmp/contigs.fa

# remove duplicates
seqkit rmdup -s tmp/contigs.fa -o tmp/contigs.fna

# cluster contigs
# cd-hit-est -i tmp/contigs.fna -o Result/02.assembly/contigs.fna -T $THREADS -M 0 -c 0.99 -d 100 -aS 0.9
mmseqs easy-cluster --split-memory-limit 0 --min-seq-id 0.95 -c 0.90 --cov-mode 0 --threads $THREADS tmp/contigs.fna  tmp/contigs tmp/ &> logs/mmseqs_contigs.log
mv tmp/contigs_rep_seq.fasta Result/02.assembly/contigs.fna

# predict
prokka Result/02.assembly/contigs.fna --force --outdir  tmp/ --prefix contigs --metagenome --cpus $THREADS --kingdom Archaea,Bacteria,Mitochondria,Viruses
prodigal -p meta -f gff -i Result/02.assembly/contigs.fna -o tmp/prodigal.gff -d tmp/prodigal.fna -a tmp/prodigal.faa 

# unique gene
# cd-hit-est -i tmp/prodigal.fna -o Result/03.gene_prediction/Unigene.fna -T 60 -M 0  -G 1 -c 0.95
mmseqs easy-cluster --split-memory-limit 0 --min-seq-id 0.95 -c 0.90 --cov-mode 0 --threads $THREADS tmp/prodigal.fna tmp/Unigene tmp/
mv tmp/Unigene_rep_seq.fasta Result/03.gene_prediction/Unigene.fna

# filter faa gtf
cut -f 1 tmp/Unigene_cluster.tsv | uniq > tmp/cluster
seqkit grep --threads $THREADS --pattern-file tmp/cluster tmp/prodigal.faa --out-file Result/03.gene_prediction/Unigene.faa
grep --fixed-strings --file tmp/cluster tmp/prodigal.gff > Result/03.gene_prediction/Unigene.gff

# count
salmon_index=tmp/salmon_index
salmon index --threads $THREADS --transcripts Result/03.gene_prediction/Unigene.fna --index $salmon_index
print_sample | parallel -j 4 salmon_quant
salmon quantmerge --quants tmp/salmon_quant/* --column TPM --genes --output tmp/salmon.tpm.tsv
salmon quantmerge --quants tmp/salmon_quant/* --column numreads --genes --output tmp/salmon.counts.tsv

# salmon2table
julia sctipts/salmon2table.jl -p $THREADS --counts -s tmp/salmon.counts.tsv.gz -o Result/03.gene_prediction/salmon.counts.xls.gz
julia sctipts/salmon2table.jl -p $THREADS -s tmp/salmon.tpm.tsv.gz -o Result/03.gene_prediction/salmon.tpm.xls.gz

# eggnog
emapper.py \
--itype proteins \
-i Result/03.gene_prediction/Unigene.faa  \
--cpu $THREADS \
-m diamond \
--data_dir /mnt/e/Alex/DATAHUB/MicroBio/eggnog \
--output_dir tmp \
--evalue 1e-03 \
--output Result/04.annotation/eggnog.annotations.xls

# GO
julia scripts/get_go_adundance.jl \
--map $MetaFlowDB/go.annotation.txt \
--eggnog Result/04.annotation/eggnog.annotations.xls.gz \
--threads 50 \
--counts Result/03.gene_prediction/salmon.counts.xls.gz \
--tpm Result/03.gene_prediction/salmon.tpm.xls.gz \
--output  /home/data/wd/Vik/X207/GO

# COG
julia scripts/get_cog_adundance.jl \
--eggnog Result/04.annotation/eggnog.annotations.xls.gz \
--counts Result/03.gene_prediction/salmon.counts.xls.gz \
--tpm Result/03.gene_prediction/salmon.tpm.xls.gz \
--output Result/04.annotation/COG_adundance \
--threads 50

tar -xf /home/data/wd/Vik/MetaFlow/MetaFlowDB/KEGG.tgz
MetaFlowDB=/home/data/wd/Vik/MetaFlow/MetaFlowDB
# KEGG
micromamba env list
micromamba activate kofamscan
exec_annotation \
--cpu 55 \
--format detail-tsv \
--profile $MetaFlowDB/KEGG/profiles/ \
--ko-list $MetaFlowDB/KEGG/ko_list \
-o Result/04.annotation/KEGG.xls \
Result/03.gene_prediction/Unigene.faa

julia scripts/get_kegg_adundance.jl \
--kofamscan Result/04.annotation/KEGG.xls \
--map $MetaFlowDB/KEGG/kegg_map.tsv \
--counts Result/03.gene_prediction/salmon.counts.xls.gz \
--tpm Result/03.gene_prediction/salmon.tpm.xls.gz \
--output Result/04.annotation/KEGG_adundance \
--threads 50

# swissprot CAZyDB TCDB VFDB CARD
source sctipts/annotation.sh

# taxonomy annotation
diamond blastp \
--db NR --taxonlist 2,2157,10239 \
--query Result/03.gene_prediction/Unigene.faa \
--threads $THREADS \
--log --evalue 1e-05 --top 10 --outfmt 100 \
--sensitive --include-lineage --out Result/06.taxonomy/NR

# diamond dda ->  megan dda
$megan/daa-meganizer \
--longReads \
--in Result/06.taxonomy/NR.daa \
--mapDB $MetaFlowDB/megan-map-Feb2022.db \
-lcaAlgorithm longReads \
--threads 55 --classify true

# megan dda -> megan rma
$megan/daa2rma \
--in Result/06.taxonomy/NR.daa \
--out Result/06.taxonomy/NR.rma \
--mapDB $MetaFlowDB/megan-map-Feb2022.db \
--longReads --lcaAlgorithm longReads \
--threads 55 --classify true \
--minScore 50 --maxExpected 0.01 --topPercent 50 

# megan rma -> table
$megan/rma2info \
--in Result/06.taxonomy/NR.rma \
--out  Result/06.taxonomy/X207.tsv \
--class2count Taxonomy --read2class Taxonomy \
--ignoreUnassigned true \
--names true --paths true --ranks true --list true --listMore true --majorRanksOnly true --verbose

# stats
seqkit fx2tab Result/02.assembly/contigs.fna -n -l -o tmp/gene_lengths.tsv
python sctipts/gene_length.py
