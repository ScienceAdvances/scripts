NJOBS=50
MEMORY=150G
mamba env create -f src/trinity.yaml
mamba activate trinity
TRINITY_HOME=$HOME/mamba/envs/trinity/opt/trinity-2.15.2

# ======== Step1 Raw sequence quality trim using fastp ========
mkidr -p tmp Resulst/00.RawData Resulst/01.CleanData Resulst/02.Assembly Resulst/03.Quantification Resulst/04.DifferentialExpression Resulst/04.Ecrihment
mkidr -p Resulst/10.Stat/fastp_report 
# ======== Step1 Raw sequence quality trim using fastp ========
ls Resulst/00.RawData/*R1.fq.gz | while read i
do
b=$(basename $i)
name=${b%%.*}
fastp -i Resulst/00.RawData/${name}.R1.fq.gz \
    -o Resulst/01.CleanData/${name}_1.fastq.gz \
    -I Resulst/00.RawData/${name}.R2.fq.gz \
    -O Resulst/01.CleanData/${name}_2.fastq.gz \
    -w 16 -h Resulst/10.Stat/fastp_report/${name}.html -j Resulst/10.Stat/fastp_report/${name}.json
done

# ======== Step2 assemble by Trinity ========
Trinity \
    --seqType fq \
    --max_memory $MEMORY  \
    --samples_file samples_Trinity.tsv \
    --CPU $NJOBS \
    --output Resulst/02.Assembly/Trinity_out

# ======== Step3 remove redundancy ========
cd-hit-est -i Resulst/02.Assembly/Trinity_out/Trinity_out.Trinity.fasta -o Resulst/02.Assembly/Trinity_out/cdhit.fa -T 30 -M 100000
$TRINITY_HOME/util/TrinityStats.pl Resulst/02.Assembly/Trinity_out/cdhit.fa > Resulst/02.Assembly/Trinity_out/TrinityStat.txt

# ======== Step3 Salmon quantification ========
$HOME/miniforge3/envs/trinity/bin/align_and_estimate_abundance.pl  \
    --transcripts Resulst/02.Assembly/Trinity_out/cdhit.fa \
    --seqType fq \
    --trinity_mode \
    --samples_file samples_Trinity.tsv \
    --est_method salmon \
    --thread_count $NJOBS \
    --prep_reference

# ======== Step4 创建表达量矩阵 ========
# 需要基因名一致（因为Step3去除了冗余）才能用用于abundance_estimates_to_matrix.pl
from Bio import SeqIO
import pandas as pd
fasta_file = "Resulst/02.Assembly/Trinity_out/cdhit.fa"
 # 遍历FASTA文件中的每个记录（序列）
rd=[]
for i in SeqIO.parse(fasta_file, "fasta"):
    rd.append(i.id)
m=pd.read_csv('Resulst/02.Assembly/Trinity_out/Trinity_out.Trinity.fasta.gene_trans_map',sep='\t',header=None)
m.loc[m[1].isin(rd),:].to_csv('Resulst/02.Assembly/Trinity_out/cdhit.fa.gene_trans_map',header=False,index=False,sep='\t')

# out: salmon.isoform.counts.matrix salmon.isoform.TPM.not_cross_norm
$HOME/miniforge3/envs/trinity/bin/abundance_estimates_to_matrix.pl \
    --est_method salmon \
    --gene_trans_map Trinity_out/cdhit.fa.gene_trans_map \
    --out_prefix Resulst/03.Quantification/salmon \
    --name_sample_by_basedir Resulst/03.Quantification/*/quant.sf

# ======== Step4 低表达量转录本的过滤 ========
$HOME/miniforge3/envs/trinity/bin/filter_low_expr_transcripts.pl \
    --matrix Resulst/03.Quantification/salmon.isoform.TPM.not_cross_norm \
    --transcripts Trinity_out/cdhit.fa \
    --min_expr_any 3 >  Resulst/02.Assembly/transcripts.fa

# ======== Step4 差异分析 ========
# 修改$TRINITY_HOME/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R中top_gene_labels_show=6
$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
    --matrix Resulst/03.Quantification/salmon_out.gene.counts.matrix \
    --method DESeq2 \
    --samples_file condition.txt \
    --contrasts compare.txt \
    --output Resulst/04.DifferentialExpression/DESeq2_results

# JOT
$TRINITY_HOME/util/misc/get_longest_isoform_seq_per_trinity_gene.pl
ls $TRINITY_HOME/Analysis/FL_reconstruction_analysis/util
$TRINITY_HOME/get_longest_isoform_seq_per_trinity_gene.pl \
../Assembly/trinity_out_dir/Trinity.fasta > longest_isoform.fasta



