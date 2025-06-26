import pandas as pd
import numpy as np
import os
import sys
import optparse


def main(meta):
    df = pd.read_csv(meta, sep='\t', header=0)
    df2 = df.loc[df.Group!="input",:]
    ncat=len(df2.Group.unique())
    if ncat < 2:
        print("Error: At least two groups are required for consensus peak analysis.")
        sys.exit(1)
    else:
        print(f"Number of groups: {ncat}")
        for x in df2.Group.unique():
            for y in ["narrowPeak","broadPeak"]:
                sample_names=df2.loc[df2.Group == x,"Sample"].values.tolist()
                ss=" ".join([f"Result/03.Peak_Calling/{y}/{x}_peaks.{y}" for x in sample_names])
                commond = f"""
                    idr --samples {ss} \
                    --input-file-type {y} \
                    --rank signal.value \
                    --output-file Result/03.Peak_Calling/{x}_peaks_{y}_IDR.bed \
                    --plot \
                    --log-output-file logs/{x}_peaks_{y}_IDR.log
                    """
                os.system(commond)

if __name__ == '__main__': 
    meta=os.environ['META']
    main(meta)
