mamba activate sc; ipython
import json
import pandas as pd
with open('genes.json') as f:
    d=json.load(f)

df=pd.read_excel('QC.xlsx')

df.columns


d={i:df[i].dropna().to_list() for i in df.columns}
with open("SCQC.json", "w") as outfile:
    json.dump(d, outfile)


from typing import Sequence,List,Union
import pandas as pd
def hsa2mmu(genelist:Sequence[str],drop:bool=False) -> Union[pd.DataFrame,List[str]]:
    from gseapy import Biomart
    q = Biomart().query(dataset='hsapiens_gene_ensembl',
        attributes=['ensembl_gene_id','external_gene_name',
                    'mmusculus_homolog_ensembl_gene',
                    'mmusculus_homolog_associated_gene_name'])
    if drop:
        return q.loc[q.external_gene_name.isin(genelist),'mmusculus_homolog_associated_gene_name'].dropna().tolist()
    else:
        return q.loc[q.external_gene_name.isin(genelist),:]

hsa2mmu(df.hsa_s.dropna(),drop=False)
hsa2mmu(df.hsa_s)



def mmu2hsa(genelist:Sequence[str],drop:bool=False):
    from gseapy import Biomart
    q = Biomart().query(dataset='mmusculus_gene_ensembl',
        attributes=['ensembl_gene_id','external_gene_name',
                    'hsapiens_homolog_ensembl_gene',
                    'hsapiens_homolog_associated_gene_name'])
    if drop:
        q.loc[q.external_gene_name.isin(genelist),'hsapiens_homolog_associated_gene_name'].dropna().tolist()
    else:
        return q.loc[q.external_gene_name.isin(genelist),:]

h2m.loc[h2m.external_gene_name.isin(['RB1']),:]
h2m.to_csv('2.csv')
hsa2mmu(['ATAD2','TP53','RB1','KRAS'],drop=False)
