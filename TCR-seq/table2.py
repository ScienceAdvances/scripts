# 统计raw fastq
micromamba activate q2
micromamba install seqkit
micromamba env list
ipython
mkdir Result/fastq
seqkit stats CleanData/*.fastq.gz --tabular --skip-err --all > Result/fastp/seqinfo_clean.tsv
import pandas as pd
import pathlib
d=pd.read_csv('Result/fastp/seqinfo_clean.tsv',sep='\t',index_col=0)
cols=['num_seqs','sum_len','Q20(%)','Q30(%)','AvgQual','GC(%)','avg_len']
d=d.loc[:,cols]

newname=['total reads','total bases','Q20(%)','Q30(%)','AvgQual','GC(%)','avg_len']
d.columns = newname
pair_end=True
d.head(2)

import numpy as np

if pair_end:
    d2=d.copy()
    d=d2.copy()
    s1=[x[0] for x in d.index.str.split(r'_\d{1,2}_L001_R[12]{1}_001.fastq.gz')]
    s2=[x.split(r'/')[-1] for x in s1]
    d2.index = np.ndarray.flatten(np.array([[f"{i}_R1",f"{i}_R2"] for i in np.unique(s2)]))
    d['Group']=s2
    d=d.groupby('Group').agg({
        'total reads':'sum',
        'total bases':'sum',
        'Q20(%)':'mean',
        'Q30(%)':'mean',
        'AvgQual':'mean',
        'GC(%)':'mean',
        'avg_len': 'mean'
        })
    d.index.names = ['Sample']
    d.loc[:,'read1 mean length'] = d2.loc[[f'{i}_R1' for i in d.index],'avg_len'].values
    d.loc[:,'read2 mean length'] = d2.loc[[f'{i}_R2' for i in d.index],'avg_len'].values
    d['Sample'] = [pathlib.Path(x[0]).stem for x in d.index.str.split('.fq')]
else:
    d['Sample'] = [pathlib.Path(x[0]).stem for x in d.index.str.split('.fa')]

d.loc[:,['Sample'] + newname].to_csv('Result/CleanFastq_Stat.xls',sep='\t',index=False)
