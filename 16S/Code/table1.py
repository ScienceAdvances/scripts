# 统计raw fastq
import pandas as pd
import pathlib
d=pd.read_csv('tmp/seqinfo.tsv',sep='\t',index_col=0)
cols=['num_seqs','sum_len','Q20(%)','Q30(%)','AvgQual','GC(%)','avg_len']
d=d.loc[:,cols]

newname=['total reads','total bases','Q20(%)','Q30(%)','AvgQual','GC(%)','avg_len']
d.columns = newname
pair_end=True
d.head(2)
if pair_end:
    d2=d.copy()
    d=d2.copy()
    s1=[x[0] for x in d.index.str.split(r'_R[12]{1}\.fastq.gz')]
    s2=[x.split(r'/')[1] for x in s1]
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
    d.loc[:,'read1 mean length'] = d2.loc[[f'RawData/{i}_R1.fastq.gz' for i in d.index],'avg_len'].values
    d.loc[:,'read2 mean length'] = d2.loc[[f'RawData/{i}_R2.fastq.gz' for i in d.index],'avg_len'].values
    d['Sample'] = [pathlib.Path(x[0]).stem for x in d.index.str.split('.fq')]
else:
    d['Sample'] = [pathlib.Path(x[0]).stem for x in d.index.str.split('.fa')]

d.loc[:,['Sample'] + newname].to_csv('Result/07.statistics/RawFastq_Stat.xls',sep='\t',index=False)
