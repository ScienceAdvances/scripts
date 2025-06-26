import pandas as pd
import numpy as np
import os
import sys
import optparse


def consensus_peak(meta):
    df = pd.read_csv(meta, sep='\t', header=0)
    df2 = df.loc[df.Group!="input",:]
    ncat=len(df2.Group.unique())
    if ncat < 2:
        print("Error: At least two groups are required for consensus peak analysis.")
        sys.exit(1)
    else:
        print(f"Number of groups: {ncat}")
        for i in df2.Group.unique():
            sample_names=df2.loc[df2.Group == i,"Sample"].values.tolist()
            bs=",".join([f"Result/03.Peak_Calling/NarrowPeak/{x}_peaks.narrowPeak" for x in sample_names[1::]])
            NarrowPeak = f"bedtools intersect -f 0.30 -r -wo -a Result/03.Peak_Calling/NarrowPeak/{sample_names[1]}_peaks.narrowPeak -b {bs} > Result/03.Peak_Calling/NarrowPeak/{i}_consensus_peaks.bed"
            os.system(NarrowPeak)
            BroadPeak = f"bedtools intersect -f 0.30 -r -wo -a Result/03.Peak_Calling/BroadPeak/{sample_names[1]}_peaks.broadPeak -b {bs} > Result/03.Peak_Calling/BroadPeak/{i}_consensus_peaks.bed"
            os.system(BroadPeak)

if __name__ == '__main__': 
    meta=os.environ['META']
    consensus_peak(meta)
