import pandas as pd
import os,re,sys

def union_genes_from_paths(paths: list):
    genes = set()
    for path in paths:
        print(path)
        temp_genes = set(get_genes_from_path(path))
        print(len(temp_genes))
        genes.update(temp_genes)
    return genes

def get_genes_from_path(path: str):
    if os.path.exists(path):
        df = pd.read_csv(path,sep='\t',index_col=0)
        return df.index.tolist()
    else:
        return []

def find_tn_pairs(readme: pd.DataFrame):
    d = dict()
    for index,row in readme.iterrows():
        c = row['Class']
        e = row['Experiment.Method']
        g = row['Group']
        if g != 'gene':
            continue
        q = row['Quant.Method']
        l = row['logTransform']
        p = row['Pathology']
        n = row['Normalized']
        path = row['Path']
        key = f"{c}_{e}_{g}_{q}_{l}_{n}"
        if key not in d:
            d[key] = dict()
        d[key][p] = path
    d2 = dict()
    for i in d:
        if len(d[i].keys()) == 2:
            d2[i] = d[i]
    return d2
