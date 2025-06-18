mamba activate sc
ipython 
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from statannotations.Annotator import Annotator
import pandas as pd
import statannotations
sns.set(style="whitegrid")

font_manager.fontManager.addfont("/usr/share/fonts/Arial/arial unicode ms.otf")
plt.rcParams['font.family'] = 'arial unicode ms'
plt.rcParams.update({'font.size': 20})

# ============================== diversity ====================================
meta=pd.read_table("Result/metadata.tsv",sep="\t",index_col=0)
df=pd.read_table("Result/postanalysis/tables.diversity.TRB.tsv",sep="\t",index_col=0)
df.index=[x[0] for x in df.index.str.split(".")]
df2 = pd.merge(meta,df,left_index=True,right_index=True)

for i in df2.columns[1::]:
    plotting_parameters = {
        'data': df2,
        'x': 'Group',
        'y': i,
        "order": ["HC",'Pre', 'Post']
    }
    fig, ax = plt.subplots(1,1,figsize=(5, 5))
    ax2=sns.boxplot(ax=ax,palette='Set2',width=0.5,**plotting_parameters)
    ax2.grid(False)
    ax2.set_xlabel("")
    ax2.set_ylabel(plotting_parameters['y'],fontsize=15)
    pairs = [
        ("HC", "Pre"),
        ("Pre", "Post"),
        ("HC", "Post")
    ]
    annot = Annotator(ax2, pairs,fontsize=10, **plotting_parameters)
    annot.configure(test='Mann-Whitney', text_format='simple', loc='inside',show_test_name=False, verbose=2)
    annot.apply_test()
    ax, test_results = annot.annotate()
    ax.tick_params(axis='both', which='major', labelsize=10)
    fig.savefig(f"Result/diversity/{i}.pdf",bbox_inches='tight')
    fig.savefig(f"Result/diversity/{i}.tiff",bbox_inches='tight')


# ============================== CDR3 Metrics ====================================
df=pd.read_table("Result/postanalysis/tables.cdr3metrics.TRB.tsv",sep="\t",index_col=0)
df.index=[x[0] for x in df.index.str.split(".")]
df2 = pd.merge(meta,df,left_index=True,right_index=True)
boxplot(df2,outdir="Result/diversity/CDR3")

def boxplot(table,outdir,group="Group",test="Mann-Whitney"):
    for i in df2.columns[1::]:
        plotting_parameters = {
            'data': table,
            'x': group,
            'y': i,
            "order": ["HC",'Pre', 'Post']
        }
        fig, ax = plt.subplots(1,1,figsize=(5, 5))
        ax2=sns.boxplot(ax=ax,palette='Set2',width=0.5,**plotting_parameters)
        ax2.grid(False)
        ax2.set_xlabel("")
        ax2.set_ylabel(plotting_parameters['y'],fontsize=15)
        pairs = [
            ("HC", "Pre"),
            ("Pre", "Post"),
            ("HC", "Post")
        ]
        annot = Annotator(ax2, pairs,fontsize=10, **plotting_parameters)
        annot.configure(test='Mann-Whitney', text_format='simple', loc='inside',show_test_name=False, verbose=2)
        annot.apply_test()
        ax, test_results = annot.annotate()
        ax.tick_params(axis='both', which='major', labelsize=10)
        fig.savefig(f"{outdir}/{i}.pdf",bbox_inches='tight')
        fig.savefig(f"{outdir}/{i}.tiff",bbox_inches='tight')

