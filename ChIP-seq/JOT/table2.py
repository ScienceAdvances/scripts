import pysam
import pandas as pd
import pathlib

def count_reads(bam_file, name):
    total_reads = 0
    mapped_reads = 0
    unique_reads = 0
    # 使用集合来记录已经出现过的read_key，加快查找判断速度
    read_positions = set()  
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile:
            total_reads += 1
            if read.is_unmapped:
                continue
            mapped_reads += 1
            read_key = (read.query_name, read.reference_name, read.reference_start)
            if read_key not in read_positions:
                read_positions.add(read_key)
                unique_reads += 1

    # 计算映射率，先以数值形式存储，后续使用时再按需格式化
    mapped_rate = mapped_reads / total_reads
    df = pd.DataFrame({
        'Sample': [name],
        'total_reads': [total_reads],
        'mapped_reads': [mapped_reads],
        'unique_reads': [unique_reads],
        'mapped_rate': [mapped_rate]
    })
    return df

names=['KO_1','KO_2','WT_2_Input','KO_2_Input','WT_2','KO_1','KO_1_Input','WT_1','WT_1_Input']

list(pathlib.Path('Result/RawBAM').glob('*_raw.bam'),names)
df = pd.DataFrame()
for x in pathlib.Path('Result/RawBAM').glob('*_raw.bam'):
    df=pd.concat([df,count_reads(x,pathlib.Path(x).stem.replace('_raw',''))])

df['mapped_rate'] = [f'{x*100:.2f}%' for x in df['mapped_rate']]
df.to_csv('Report/table2.csv',index=False)
