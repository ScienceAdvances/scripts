import tempfile
import os
import pathlib
import shutil
from pathlib import Path
import json
import numpy as np
import pandas as pd
import biom
import qiime2
from qiime2 import Metadata
from q2_types.per_sample_sequences import PairedEndFastqManifestPhred33V2
import qiime2.plugins.dada2.actions as dada2_actions
import qiime2.plugins.feature_table.actions as feature_table_actions
import qiime2.plugins.metadata.actions as metadata_actions
import qiime2.plugins.feature_classifier.actions as feature_classifier_actions
import qiime2.plugins.taxa.actions as taxa_actions
from qiime2.plugins.cutadapt.methods import trim_paired
import qiime2.plugins.phylogeny.actions as phylogeny_actions
from qiime2.plugins.vsearch.methods import cluster_features_de_novo
import Bio
import os

THREADS = os.environ['THREADS']
MiroFlowDB = os.environ['MiroFlowDB']
META = os.environ['META']
SEQ="manifest.tsv"
CLASSIFIER=f'{MiroFlowDB}/silva-138-99-nb-classifier.qza'

# ======================= import data =======================
sample_metadata_md = Metadata.load(META)
classifier=qiime2.Artifact.load(CLASSIFIER)
demultiplexed_sequences = qiime2.Artifact.import_data(
    'SampleData[PairedEndSequencesWithQuality]',
    SEQ,
    PairedEndFastqManifestPhred33V2,
)
# ======================= 2. 生成特征表和代表序列 =======================
feature_table, asv_sequences, dada2_stats = dada2_actions.denoise_paired(
    demultiplexed_seqs=demultiplexed_sequences,
    trunc_len_f=0,
    trunc_len_r=0,
    trim_left_f=0,
    trim_left_r=0,
    max_ee_f=1.2,
    max_ee_r=1.2,
    min_overlap=12,
    pooling_method='independent',
    chimera_method='consensus',
    n_threads=THREADS,
)

feature_table.save("Result/02.asv_taxonomy/feature_table.qza")
asv_sequences.save("Result/02.asv_taxonomy/RepresentSequence.qza")

# Reviewing the DADA2 run statistics
dada2_stats.view(Metadata).to_dataframe().to_csv("Result/02.asv_taxonomy/DADA2_stats.xls",sep='\t')

# 3. feature table and feature data
feature_table1, = feature_table_actions.filter_samples(
    table=feature_table,
    min_frequency=0,
    min_features=0,
    metadata=sample_metadata_md
)

filtered_sequences_1, = feature_table_actions.filter_seqs(
    data=asv_sequences,
    table=feature_table,
)
# classify
taxonomy, = feature_classifier_actions.classify_sklearn(
    classifier=classifier,
    reads=filtered_sequences_1,
)

# ASVs with an absolute count sum less than 50 and prevalence below 1% were excluded
# 1 Filtering samples with low sequence counts 过滤稀有ASV
filtered_table_1, = feature_table_actions.filter_features(
    table=feature_table1,
    min_frequency=50,
    min_samples=1
)
# 2 过滤线粒体和叶绿体
filtered_table_2, = taxa_actions.filter_table(
    table=filtered_table_1,
    taxonomy=taxonomy,
    mode='contains',
    include='p__',
    exclude='p__;,Chloroplast,Mitochondria',
)

filtered_table_2.view(biom.Table).to_dataframe().shape
# 4）过滤有较少的ASV总量的sample
# filtered_table_3, = feature_table_actions.filter_samples(
#     table=filtered_table_2,
#     min_frequency=10000,
# )
filtered_table_2.view(biom.Table).to_dataframe().shape

# 更新特征序列
filtered_sequences_2, = feature_table_actions.filter_seqs(
    data=filtered_sequences_1,
    table=filtered_table_2,
)
filtered_sequences_2.export_data("MiroFlow/tmp")

# ======================= ASV ID 重命名 =======================
ftable = filtered_table_2.view(biom.Table).to_dataframe()
idmaps = ftable.copy()
idmaps['ASVID'] = [f'ASV{x}' for x in range(1,ftable.shape[0]+1)]
idmaps = idmaps.loc[:,['ASVID']]
table=pd.merge(idmaps,ftable,left_index=True,right_index=True)
table.to_csv("Result/02.asv_taxonomy/ASV.xls",sep='\t',index=False)
table=pd.read_csv("Result/02.asv_taxonomy/ASV.xls",sep='\t')
idmaps.to_csv("tmp/idmaps.xls",sep='\t',index=True)

# ======================= 处理特征序列 =======================
os.system("julia Code/src/rename_fa.jl")

# ======== 处理特征表 ============
ftable = qiime2.Artifact.import_data(
    type='FeatureTable[Frequency]',
    view=table.set_index('ASVID').T
)

seq = qiime2.Artifact.import_data(
    'FeatureData[Sequence]',
    'Result/02.asv_taxonomy/RepresentSequence.fasta'
)
# 重新注释
taxonomy_2, = feature_classifier_actions.classify_sklearn(
    classifier=classifier,
    reads=seq,
)
tax_df=taxonomy_2.view(pd.DataFrame)
tax_text = tax_df.Taxon.str.split(';')
columns=["Kingdom","Phylum","Class","Order","Family","Genus","Sepcies"]
tax=pd.DataFrame(columns=columns)
for i in tax_text:
    i.extend(["__"]*(7-len(i)))
    t1 = pd.DataFrame(i,index=columns).replace("__","").T
    tax=pd.concat([tax, t1], axis=0)
tax.index = tax_df.index
tax.index.names = ["ASVID"]
tax['Confidence'] = tax_df.Confidence
tax.to_csv("Result/02.asv_taxonomy/silva_taxonomy.xls",sep='\t')

for i in range(1,8):
    a, = taxa_actions.collapse(
        table=ftable,
        taxonomy=taxonomy_2,
        level=i,
    )
    a2=a.view(biom.Table).to_dataframe()
    a2.index.names = ['Tax']
    a2.to_csv(f"Result/02.asv_taxonomy/silva_level_{i}.xls",sep='\t')


# Phylogenetic tree construction
phylogeny = phylogeny_actions.align_to_tree_mafft_fasttree(
    sequences=seq,
    n_threads=THREADS,
)

with tempfile.TemporaryDirectory() as tmp_dir:
    phylogeny.tree.export_data(tmp_dir)
    shutil.copy2(pathlib.Path(tmp_dir).joinpath('tree.nwk'), pathlib.Path('Result/02.asv_taxonomy').joinpath('unrooted_tree.nwk'))

with tempfile.TemporaryDirectory() as tmp_dir:
    phylogeny.rooted_tree.export_data(tmp_dir)
    shutil.copy2(pathlib.Path(tmp_dir).joinpath('tree.nwk'), pathlib.Path('Result/02.asv_taxonomy').joinpath('rooted_tree.nwk'))
