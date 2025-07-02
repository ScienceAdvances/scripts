#!/usr/bin/bash
# ========================== mkdir ===========================
mkdir -p tmp logs Result/01.clean_data  Result/02.asv_taxonomy  Result/03.alpha_diversity Result/04.beta_diversity Result/05.differential_abundance Result/06.functional_prediction Result/07.statistics/fastp_report


# ========================== global config ===========================
export WORKDIR=/home/zktbk0/Alex/MiroFlow
export MiroFlowDB=/home/zktbk0/Alex/MiroFlow/MiroFlowDB
export META="meta.tsv"
export THREADS=32
mamba activate MiroFlow
cd $WORKDIR
source Code/src/fcuntions.sh

# ========================== install conda env ===========================
# 安装qiime2x https://library.qiime2.org/quickstart/amplicon
mamba env create -f Code/src/MiroFlow.yaml
mamba activate MiroFlow
# 或者使用容器
apptainer build --sandbox amplicon docker://quay.io/qiime2/amplicon:2025.04

# ========================== install R packages ===========================
Rscript -e 'pkgs=list("microeco","ChiLiubio/mecodev","ScienceAdvances/using","ricardo-bion/ggradar","gmteunisse/ggnested");lapply(pkgs,BiocManager::install)'
# remotes::install_local("./ggpicrust2-main")
# BiocManager::install("cafferychen777/ggpicrust2")

# ================== get samples ================
mapfile -t samples < <(awk 'NR>1 {print $1}' $META)
export samples

# ================= clean fastq files ================
export -f run_fastp
print_sample | parallel -j 16 run_fastp
seqkit stats Result/01.clean_data/*.fq.gz --tabular --skip-err --all > tmp/seqinfo_clean.tsv
seqkit stats Result/00.raw_data/*.fq.gz --tabular --skip-err --all > tmp/seqinfo_raw.tsv

# ================== run QIIME2 ================
python Code/01.qiime2.py

# ================== run userach ================
# Annotation
usearch11=$MiroFlowDB/usearch11.0.667_i86linux64
usearch10=$MiroFlowDB/usearch10.0.240_i86linux64

$usearch11 \
-sintax Result/02.asv_taxonomy/RepresentSequence.fasta \
-sintax_cutoff 0.8 \
-db $MiroFlowDB/rdp_16s_v18.fa \
-strand both \
-tabbedout Result/02.asv_taxonomy/rdp_taxonomy.xls \
-threads $THREADS

# filter taxonomy
julia Code/src/tax_filter.jl

#usearch analysis
$usearch11 -otutab_stats Result/02.asv_taxonomy/ASV.xls -output  Result/02.asv_taxonomy/ASV_use_report.xls
$usearch11 -calc_distmx Result/02.asv_taxonomy/RepresentSequence.fasta -tabbedout Result/02.asv_taxonomy/ASV_distmx.xls

#core ASVs
$usearch11 -otutab_core Result/02.asv_taxonomy/ASV.xls -distmxin Result/02.asv_taxonomy/ASV_distmx.xls -tabbedout Result/02.asv_taxonomy/core_ASV.xls

#alpha diversity
$usearch11 -otutab_rare Result/02.asv_taxonomy/ASV.xls -sample_size 10000 -output Result/02.asv_taxonomy/ASVtab_uparse_10k.xls
$usearch11 -alpha_div Result/02.asv_taxonomy/ASVtab_uparse_10k.xls -output Result/03.alpha_diversity/alpha_ASV_10k.xls
$usearch10 -alpha_div_rare Result/02.asv_taxonomy/ASVtab_uparse_10k.xls -output Result/03.alpha_diversity/alpha_rare.xls

#beta_diversity
$usearch11 -cluster_aggd  Result/02.asv_taxonomy/ASV_distmx.xls -treeout Result/04.beta_diversity/ASV_uparse.tree
$usearch11 -beta_div Result/03.alpha_diversity/ASVtab_uparse_10k.xls -tree Result/04.beta_diversity/ASV_uparse.tree -filename_prefix Result/04.beta_diversity/ -metrics bray_curtis,jaccard,unifrac,unifrac_binary


# ================== run picrust2 ================
mamba activate picrust2
picrust2_pipeline.py \
    --study_fasta Result/02.asv_taxonomy/RepresentSequence.fasta \
    --input Result/02.asv_taxonomy/ASV.xls \
    --output Result/06.functional_prediction \
    --min_align 0.8 \
    --processes $THREADS

