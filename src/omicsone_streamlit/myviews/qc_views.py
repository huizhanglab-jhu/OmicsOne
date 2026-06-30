import numpy as np
import pandas as pd
import streamlit as st
import re,sys,os
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from tqdm import tqdm


from dataclasses import dataclass, field
from typing import List, Dict, Literal, Iterable, Optional, Union


Pathology = Literal["Tumor", "NAT", "Mixed"]
@dataclass
class QCInput:
    omics_select: str
    data: pd.DataFrame
    meta: pd.DataFrame
    sel_set: str
    sel_group: str
    pathology_type: Pathology = "Tumor"
    
    


def corrfunc(x,y,**kws):
    r, _ = stats.spearmanr(x,y,nan_policy='omit')
    ax = plt.gca()
    ax.annotate("r = {:.2f}".format(r), xy=(.1,.5),
                xycoords = ax.transAxes, fontsize=30)

def qc_corr_page(samples, data):

    cols = st.columns([1,1])
    with cols[1]:
        sample_selects = st.multiselect("Samples",samples)

        data2 = data.loc[:, sample_selects]
        data2 = data2.replace(0,np.nan)
        # data2 = data2.dropna()
        st.write('features: {:d}'.format(data2.shape[0]))
        

    with cols[0]:
        if len(sample_selects) >= 2:
            g = sns.PairGrid(data2, palette=["red"])
            g.map_lower(sns.scatterplot, s=10, legend=False)
            g.map_diag(sns.histplot, kde=False)
            g.map_upper(corrfunc)

            st.pyplot(g.fig)
            
def get_glycoform(index):
    # glycoform, sequence, protein, site = index
    glycoform, proteins = index
    return glycoform

def get_glyco_gene(index):
    # glycoform, sequence, protein, site = index
    glycoform, proteins = index
    protein_list = proteins.split(";")
    items = protein_list[0].split("|")
    return items[1]

def get_glyco_seq(index):
    # glycoform, sequence, protein, site = index
    glycoform, proteins = index
    sequence = glycoform.split("-")[0]
    return sequence

def get_glycosite(index):
    # glycoform, sequence, protein, site = index
    glycoform, proteins = index
    protein_list = proteins.split(";")
    items = protein_list[0].split("|")
    return f"{items[0]}@{items[4]}"

def get_phospho_site(phospho_idx):
    items = phospho_idx[0].split("|")
    protein = items[1]
    site = phospho_idx[0].split("_")[-1]
    return f"{protein}@{site}" 

def get_phospho_gene(phospho_idx):
    gene = phospho_idx[1]
    return gene

def get_phospho_seq(phospho_idx):
    seq = phospho_idx[2]
    return seq

def get_protein_gene(protein_idx):
    items = protein_idx.split("|")
    gene = items[0]
    return gene

def get_protein_id(protein_idx):
    items = protein_idx.split("|")
    if len(items) == 2:
        protein = items[1]
    else:
        protein = items[0]
    return protein


@st.cache_data
def build_uniq_feature(df, meta_df, sel_set,sel_feature="phospho_seq" ):
    # read sample file
    samples = list(df.columns.values)
    
    meta_df = meta_df[meta_df.index.isin(samples)]
    
    uniq_sets = meta_df[sel_set].unique()
    
    sample_set_map = dict([(index,row[sel_set]) for index,row in meta_df.iterrows()])

    df = df.replace('None',np.nan) \
            .map(float) \
            .dropna(how='all')
    
    set_feature_map = dict([(i,set()) for i in uniq_sets])

    for sample in tqdm(sample_set_map):
        idx_list = list(df[pd.notna(df[sample])].index)
        features = []
        if sel_feature == "phospho_seq":
            features = [get_phospho_seq(idx) for idx in idx_list]
        elif sel_feature == "phospho_protein":
            features = [get_phospho_gene(idx) for idx in idx_list]
        elif sel_feature == "protein_gene":
            features = [get_protein_gene(idx) for idx in idx_list]
        elif sel_feature == "protein_id":
            features = [get_protein_id(idx) for idx in idx_list]
        elif sel_feature == "glycoform":
            features = [get_glycoform(idx) for idx in idx_list]
        sample_set = sample_set_map[sample]
        set_feature_map[sample_set].update(set(features))
    
        
    rows = []
    pep50 = []
    set_size = len(set_feature_map)
    for sample_set in tqdm(sorted(set_feature_map)):
        pep = set_feature_map[sample_set]
        uniqs = 0
        
        common_50 = 0
        common_100 = 0
        not_uniq = 0
        
        for p in pep:
            shares = 0
            uniq = True
            for s in sorted(set_feature_map):
                if s != sample_set and p in set_feature_map[s]:
                    shares += 1
                    uniq = False
            if uniq:
                uniqs += 1 
            else:
                if shares + 1 == set_size:
                    common_100 += 1
                elif shares + 1 >= set_size/2:
                    common_50 += 1
                    pep50.append(p)
                not_uniq += 1
            # print(p,shares,sample_size)
        # break
        rows.append([sample_set, uniqs,not_uniq, common_50, common_100])

    result = pd.DataFrame(rows, columns=['SampleSet', 'uniq', 'not_uniq','common_50','common_100'])

    return result
    
def statistics_view(cfg: QCInput):
    omics_select = cfg.omics_select
    if omics_select == 'Phosphoproteomics':
        phospho_statistics_view(cfg)
    elif omics_select == 'Proteomics':
        proteomics_statistics_view(cfg)
    elif omics_select == 'Glycoproteomics':
        glycoproteomics_statistics_view(cfg)

def proteomics_statistics_view(cfg: QCInput):
    df = cfg.data
    meta_df = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    left_col, right_col = st.columns(2)
    
    with left_col:  
        num_rows = df.shape[0]
        num_cols = df.shape[1]
        
        num_genes = len(set([get_protein_gene(idx) for idx in df.index]))
        num_proteins = len(set([get_protein_id(idx) for idx in df.index]))

        # Create a DataFrame for the statistics
        headers = ['Table.Rows', 'Table.Columns', 'ID.Genes', 'ID.Proteins']
        values  = [num_rows, num_cols, num_genes, num_proteins]
        
        groups = meta_df[sel_group].unique()  
        for g in sorted(groups):
            samples = list(meta_df[meta_df[sel_group] == g].index)
            samples = [i for i in samples if i in df.columns]
            headers.append(f'Sample.{g}')
            values.append(len(samples))

        stats_data = {
            'Item': headers,
            'Value': values
        }
        stats_df = pd.DataFrame(stats_data)
        
        st.dataframe(stats_df)
        
def glycoproteomics_statistics_view(cfg: QCInput):
    df = cfg.data
    meta_df = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    # st.write("Glycoproteomics Statistics")
    left_col, right_col = st.columns(2)
    
    with left_col:
        num_rows = df.shape[0]
        num_cols = df.shape[1]
        
        num_glycoforms = len(set([get_glycoform(idx) for idx in df.index]))
        num_glycosites = len(set([get_glycosite(idx) for idx in df.index]))
        num_genes = len(set([get_glyco_gene(idx) for idx in df.index]))
        num_seqs = len(set([get_glyco_seq(idx) for idx in df.index]))

        # Create a DataFrame for the statistics
        headers = ['Table.Rows', 'Table.Columns','ID.Glycoforms', 
                   'ID.Glycosites', 'ID.Genes', 'ID.Sequences']
        values  = [num_rows, num_cols, num_glycoforms, 
                   num_glycosites, num_genes, num_seqs]
        
        groups = meta_df[sel_group].unique()  
        # st.write(groups)
        # st.write(sel_group)
        for g in sorted(groups):
            # st.write(g)
            samples = list(meta_df[meta_df[sel_group] == g].index)
            # if pathology_type == "Mixed":
            #     suffix_map = {"Tumor": "_T", "NAT": "_N", "Replicate": "_R"}
            #     suffix = suffix_map.get(g)
            #     if suffix:  # 只在 g 合法时改名
            #         samples = [f"{s}{suffix}" for s in samples]
            samples = [i for i in samples if i in df.columns]
            headers.append(f'Sample.{g}')
            values.append(len(samples))

        stats_data = {
            'Item': headers,
            'Value': values
        }
        stats_df = pd.DataFrame(stats_data)
        
        st.dataframe(stats_df)
    
def phospho_statistics_view(cfg: QCInput):
    df = cfg.data
    meta_df = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    left_col, right_col = st.columns(2)
    
    with left_col:
        num_rows = df.shape[0]
        num_cols = df.shape[1]
        
        num_phosphosites = len(set([get_phospho_site(idx) for idx in df.index]))
        num_genes = len(set([get_phospho_gene(idx) for idx in df.index]))
        num_seqs = len(set([get_phospho_seq(idx) for idx in df.index]))

        # Create a DataFrame for the statistics
        headers = ['Table.Rows', 'Table.Columns', 'ID.Phosphosites', 'ID.Genes', 'ID.Sequences']
        values  = [num_rows, num_cols, num_phosphosites, num_genes, num_seqs]
        
        groups = meta_df[sel_group].unique()  
        for g in sorted(groups):
            samples = list(meta_df[meta_df[sel_group] == g].index)
            samples = [i for i in samples if i in df.columns]
            headers.append(f'Sample.{g}')
            values.append(len(samples))

        stats_data = {
            'Item': headers,
            'Value': values
        }
        stats_df = pd.DataFrame(stats_data)
        
        st.dataframe(stats_df)
        
def reproducibility_view(cfg: QCInput):
    omics_select = cfg.omics_select
    if omics_select == 'Phosphoproteomics':
        phosphoproteomics_reproducibility_view(cfg)
    elif omics_select == 'Proteomics':
        proteomics_reproducibility_view(cfg)
    elif omics_select == 'Glycoproteomics':
        glycoproteomics_reproducibility_view(cfg)

def phosphoproteomics_reproducibility_view(cfg: QCInput):
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set

    st.write("Phosphoproteomics Reproducibility")
    feature_type = "phospho_seq"
    result = build_uniq_feature(data.copy(deep=True), meta.copy(deep=True), sel_set,sel_feature=feature_type)
    with st.expander("Data"):
        st.dataframe(result)
    with st.expander("Settings"):
        figsize_cols = st.columns([1,1,4])
        with figsize_cols[0]:
            fig_width = st.number_input("Figure.Width", min_value=1, value=10)
        with figsize_cols[1]:
            fig_height = st.number_input("Figure.Height", min_value=1, value=8)
    from plots.bars import plot_uniq
    fig = plot_uniq(result, feature_type = feature_type, figsize=(fig_width,fig_height), title = "Phosphoproteomics Reproducibility")    
    st.pyplot(fig)
    
def proteomics_reproducibility_view(cfg: QCInput):
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set
    st.write("Proteomics Reproducibility")
    feature_type = "protein_gene"
    result = build_uniq_feature(data.copy(deep=True), meta.copy(deep=True), sel_set, sel_feature=feature_type)
    with st.expander("Data"):
        st.dataframe(result)
    with st.expander("Settings"):
        figsize_cols = st.columns([1,1,4])
        with figsize_cols[0]:
            fig_width = st.number_input("Figure.Width", min_value=1, value=10)
        with figsize_cols[1]:
            fig_height = st.number_input("Figure.Height", min_value=1, value=8)
    from plots.bars import plot_uniq
    fig = plot_uniq(result, feature_type=feature_type, figsize=(fig_width,fig_height), title = "Proteomics Reproducibility")    
    st.pyplot(fig)
    
def glycoproteomics_reproducibility_view(cfg: QCInput):
    st.write("Glycoproteomics Reproducibility")
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set
    feature_type = "glycoform"
    result = build_uniq_feature(data.copy(deep=True), 
                                meta.copy(deep=True), sel_set, 
                                sel_feature=feature_type)
    with st.expander("Data"):
        st.dataframe(result)
    with st.expander("Settings"):
        figsize_cols = st.columns([1,1,4])
        with figsize_cols[0]:
            fig_width = st.number_input("Figure.Width", min_value=1, value=10)
        with figsize_cols[1]:
            fig_height = st.number_input("Figure.Height", min_value=1, value=8)
    from plots.bars import plot_uniq
    fig = plot_uniq(result, feature_type=feature_type, figsize=(fig_width,fig_height), title = "Glycoproteomics Reproducibility")    
    st.pyplot(fig)
    
def pca_view(cfg: QCInput):
    omics_select = cfg.omics_select
    if omics_select == 'Phosphoproteomics':
        phosphoproteomics_pca_view(cfg)
    elif omics_select == 'Proteomics':
        proteomics_pca_view(cfg)
    elif omics_select == 'Glycoproteomics':
        glycoproteomics_pca_view(cfg)
        
def phosphoproteomics_pca_view(cfg: QCInput):
    st.write("Phosphoproteomics PCA")
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    left_col, right_col = st.columns(2)
    sample_group_map = dict(zip(meta.index, meta[sel_group]))
    cols = [i for i in data.columns.values if i != "Intensity.Reference"]
    data_df = data[cols]
    
    data_df =  data_df.replace([np.inf, -np.inf],np.nan)
    data_df = data_df.replace('None',np.nan) \
        .map(float) \
        .dropna()
    data_df = data_df.T
    
    if isinstance(data_df.columns, pd.MultiIndex):
        data_df.columns = data_df.columns.map(lambda t: "@".join(map(str, t)))
    
    data_df[sel_group] = [sample_group_map[index] for index,row in data_df.iterrows()]

    from plots.pca import plot_pca
    with left_col:
        with st.spinner('Running... Please wait...'):
            fig, df_pca = plot_pca(data_df)
            st.dataframe(df_pca)
    with right_col:
        st.pyplot(fig)

def proteomics_pca_view(cfg: QCInput):
    st.write("Proteomics PCA")
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    left_col, right_col = st.columns(2)
    sample_group_map = dict(zip(meta.index, meta[sel_group]))
    cols = [i for i in data.columns.values if i != "Intensity.Reference"]
    data_df = data[cols]
    data_df = data_df.replace('None',np.nan) \
        .map(float) \
        .dropna()
    data_df = data_df.T
    
    if isinstance(data_df.columns, pd.MultiIndex):
        data_df.columns = data_df.columns.map(lambda t: "@".join(map(str, t)))
        
    data_df[sel_group] = [sample_group_map[index] for index,row in data_df.iterrows()]

    from plots.pca import plot_pca
    with left_col:
        with st.spinner('Running... Please wait...'):
            fig,df_pca = plot_pca(data_df)
        st.dataframe(df_pca)
    with right_col:
        st.pyplot(fig)
    
def glycoproteomics_pca_view(cfg: QCInput):
    st.write("Glycoproteomics PCA")
    data = cfg.data
    meta = cfg.meta
    sel_set = cfg.sel_set
    sel_group = cfg.sel_group
    left_col, right_col = st.columns(2)
    sample_group_map = dict(zip(meta.index, meta[sel_group]))
    cols = [i for i in data.columns.values if i != "Intensity.Reference"]
    data_df = data[cols]
    data_df = data_df.replace('None',np.nan) \
        .map(float) \
        .dropna()
    data_df = data_df.T
    # flat the columns names if multi-index
    if isinstance(data_df.columns, pd.MultiIndex):
        data_df.columns = data_df.columns.map(lambda t: "@".join(map(str, t)))
    data_df[sel_group] = [sample_group_map[index] for index,row in data_df.iterrows()]

    from plots.pca import plot_pca
    with left_col:
        with st.spinner('Running... Please wait...'):
            fig,df_pca = plot_pca(data_df)
        st.dataframe(df_pca)
    with right_col:
        st.pyplot(fig)
    
    
# missing values
def missing_values_view(cfg: QCInput):
    omics_select = cfg.omics_select
    if omics_select == 'Phosphoproteomics':
        phosphoproteomics_missing_values_view(cfg)
    elif omics_select == 'Proteomics':
        proteomics_missing_values_view(cfg)
    elif omics_select == 'Glycoproteomics':
        glycoproteomics_missing_values_view(cfg)
        
def phosphoproteomics_missing_values_view(cfg: QCInput):
    st.write("Phosphoproteomics Missing Values")
    
def proteomics_missing_values_view(cfg: QCInput):
    st.write("Proteomics Missing Values")
    
def glycoproteomics_missing_values_view(cfg: QCInput):
    st.write("Glycoproteomics Missing Values")