/usr/bin/bash
wget https://www.encodeproject.org/files/ENCFF647ULQ/@@download/ENCFF647ULQ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF111LIY/@@download/ENCFF111LIY.fastq.gz
wget https://www.encodeproject.org/files/ENCFF319MIE/@@download/ENCFF319MIE.fastq.gz
wget https://www.encodeproject.org/files/ENCFF217LJP/@@download/ENCFF217LJP.fastq.gz

# ============================ racoon_clip ================================
wget https://github.com/ZarnackGroup/racoon_clip/archive/refs/tags/v1.1.3.zip
unzip v1.1.3.zip
conda install -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*'
mamba create -n racoon_clip python=3.9.0 pip
mamba install mamba -c conda-forge
mamba activate racoon_clip
cd racoon_clip-1.1.3
pip install -e .
pip install pulp==2.7.0
cd ..
racoon_clip run --configfile configfile.yaml --cores 32

# ============================ iCount ================================
mkdir Result/iCount
mamba activate racoon_clip
mamba install icount -c bioconda
iCount -h
for i in ENCFF111LIY ENCFF217LJP ENCFF319MIE ENCFF647ULQ
do
iCount xlsites \
Result/aligned/${i}.Aligned.sortedByCoord.out.bam \
Result/iCount/${i}.cDNA_unique.bed  \
Result/iCount/${i}.cDNA_multiple.bed \
Result/iCount/${i}.cDNA_skipped.bam \
--group_by start \
--quant cDNA
done

# ============================ clippy ================================
mkdir Result/clippy
mamba env create -f script/clippy.yaml
mamba activate clippy
for i in ENCFF111LIY ENCFF217LJP ENCFF319MIE ENCFF647ULQ
do
clippy \
--thread 32 \
--window_size 10 \
--intergenic_peak_threshold 5 \
--input_bed Result/iCount/${i}.cDNA_unique.bed \
--output_prefix Result/clippy/${i} \
--annotation gencode.v46.chr_patch_hapl_scaff.annotation.gtf \
--genome_file GRCh38.p14.genome.fa.fai
done

# ============================ PEKA ================================
mamba env create -f script/PEKA.yaml
mamba activate peka
# Install the latest version of PEKA from the main branch
pip install git+https://github.com/ulelab/peka@main
pip install webcolors
export TMPDIR='Result/tmp'
mkdir Result/PEKA

for i in ENCFF111LIY ENCFF217LJP ENCFF319MIE ENCFF647ULQ
do
cp Result/clippy/${i}_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed Result/PEKA/${i}.peak.bed
python script/peka.py \
--inputpeaks Result/PEKA/${i}.peak.bed \
--inputxlsites Result/iCount/${i}.cDNA_unique.bed \
--genomefasta GRCh38.p14.genome.fa \
--genomeindex GRCh38.p14.genome.fa.fai \
--regions gencode.v46.chr_patch_hapl_scaff.annotation.gtf \
--outputpath Result/PEKA \
--specificregion other_exon \
--repeats unmasked \
--kmerlength 5 \
--smoothing 6 \
--percentile 0.7 \
--topn 20 \
--distalwindow 150 \
--window 25
cut -f 1,9-110 Result/PEKA/${i}.cDNA_unique_5mer_distribution_other_exon.tsv > Result/PEKA/${i}.rtxn.tsv
python script/plotRelativeOccurrenceHeatmap.py \
    Result/PEKA/${i}.cDNA_unique_5mer_distribution_other_exon.tsv \
    Result/PEKA/${i}.rtxn.tsv \
    Result/PEKA/ \
    20 \
    40
done

# 因为位点太少 富集不到 所以下调阈值percentile为0.3
i='ENCFF647ULQ'
cp Result/clippy/${i}_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed Result/PEKA/${i}.peak.bed
python script/peka.py \
--inputpeaks Result/PEKA/${i}.peak.bed \
--inputxlsites Result/iCount/${i}.cDNA_unique.bed \
--genomefasta GRCh38.p14.genome.fa \
--genomeindex GRCh38.p14.genome.fa.fai \
--regions gencode.v46.chr_patch_hapl_scaff.annotation.gtf \
--outputpath Result/PEKA \
--specificregion other_exon \
--repeats unmasked \
--kmerlength 5 \
--smoothing 6 \
--percentile 0.3 \
--topn 20 \
--distalwindow 150 \
--window 25
cut -f 1,9-110 Result/PEKA/${i}.cDNA_unique_5mer_distribution_other_exon.tsv > Result/PEKA/${i}.rtxn.tsv
python script/plotRelativeOccurrenceHeatmap.py \
    Result/PEKA/${i}.cDNA_unique_5mer_distribution_other_exon.tsv \
    Result/PEKA/${i}.rtxn.tsv \
    Result/PEKA/ \
    20 \
    40