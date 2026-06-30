import streamlit as st
from utils.fasta import get_gene_map
import os, re,sys
import pandas as pd
import numpy as np
from utils.omicsone_pca import pca_view
from utils.omicsone_umap import umap_view
from utils.omicsone_tsne import tsne_view
from pathlib import Path


@st.cache_data
def perform_get_gene_map(fasta_path):
    gene_map = get_gene_map(fasta_path)
    return gene_map


def single_pca(gene_map, readme, data_dir, job_dir):

    file_cols = st.columns([1,2])
    with file_cols[0]:
        omics_list = [i for i in readme['Class'].unique() if not re.search('^Other',i)]
        omics_select = st.selectbox('omics', omics_list)

    with file_cols[1]:
        file_list = readme[readme['Class']==omics_select]['Path'].unique()
        file_select = st.selectbox('file', file_list)

    meta_cols = st.columns([1,2])
    with meta_cols[0]:
        meta_path_list = readme[readme['Class'].str.contains('Other_Meta')]['Path'].unique()
        meta_path_list = [i for i in meta_path_list if not re.search('CaseList',i)]
        meta_path_select = st.selectbox('meta_path', meta_path_list, key='meta_path')

    path_file = os.path.join(data_dir,file_select)
    path_meta = os.path.join(data_dir,meta_path_select)
    data_df = pd.read_csv(path_file, sep='\t',index_col=0)
    meta_df = pd.read_csv(path_meta, sep='\t',index_col=0, header=[0,1])
    meta_df = meta_df.replace(np.nan, 'NA')
    
    
    with meta_cols[1]:
        meta_options = meta_df.columns.tolist()
        meta_options = [i for i in meta_options if i[1] in ['ORD','BIN','NOM']]
        meta_select = st.selectbox('meta', meta_options, key='meta')

    with st.expander("Data"):
        tabs = st.tabs(["Data", "Meta"])
        with tabs[0]:
            st.dataframe(data_df)
        with tabs[1]:
            st.dataframe(meta_df)
    
    with st.expander("Figure"):
        pca_view(data_df, meta_df, meta_select, job_dir)


def single_umap(gene_map, readme, data_dir, job_dir):

    file_cols = st.columns([1,2])
    with file_cols[0]:
        omics_list = [i for i in readme['Class'].unique() if not re.search('^Other',i)]
        omics_select = st.selectbox('omics', omics_list, key='umap_omics')

    with file_cols[1]:
        file_list = readme[readme['Class']==omics_select]['Path'].unique()
        file_select = st.selectbox('file', file_list, key='umap_file')

    meta_cols = st.columns([1,2])
    with meta_cols[0]:
        meta_path_list = readme[readme['Class'].str.contains('Other_Meta')]['Path'].unique()
        meta_path_list = [i for i in meta_path_list if not re.search('CaseList',i)]
        meta_path_select = st.selectbox('meta_path', meta_path_list, key='umap_meta_path')

    path_file = os.path.join(data_dir,file_select)
    path_meta = os.path.join(data_dir,meta_path_select)
    data_df = pd.read_csv(path_file, sep='\t',index_col=0)
    meta_df = pd.read_csv(path_meta, sep='\t',index_col=0, header=[0,1])
    meta_df = meta_df.replace(np.nan, 'NA')
    
    
    with meta_cols[1]:
        meta_options = meta_df.columns.tolist()
        meta_options = [i for i in meta_options if i[1] in ['ORD','BIN','NOM']]
        meta_select = st.selectbox('meta', meta_options, key='umap_meta')

    with st.expander("Data"):
        tabs = st.tabs(["Data", "Meta"])
        with tabs[0]:
            st.dataframe(data_df)
        with tabs[1]:
            st.dataframe(meta_df)
    
    with st.expander("Figure"):
        umap_view(data_df, meta_df, meta_select, job_dir)



def single_tsne(gene_map, readme, data_dir, job_dir):

    file_cols = st.columns([1,2])
    with file_cols[0]:
        omics_list = [i for i in readme['Class'].unique() if not re.search('^Other',i)]
        omics_select = st.selectbox('omics', omics_list, key='tsne_omics')

    with file_cols[1]:
        file_list = readme[readme['Class']==omics_select]['Path'].unique()
        file_select = st.selectbox('file', file_list, key='tsne_file')

    meta_cols = st.columns([1,2])
    with meta_cols[0]:
        meta_path_list = readme[readme['Class'].str.contains('Other_Meta')]['Path'].unique()
        meta_path_list = [i for i in meta_path_list if not re.search('CaseList',i)]
        meta_path_select = st.selectbox('meta_path', meta_path_list, key='tsne_meta_path')

    path_file = os.path.join(data_dir,file_select)
    path_meta = os.path.join(data_dir,meta_path_select)
    data_df = pd.read_csv(path_file, sep='\t',index_col=0)
    meta_df = pd.read_csv(path_meta, sep='\t',index_col=0, header=[0,1])
    meta_df = meta_df.replace(np.nan, 'NA')
    
    
    with meta_cols[1]:
        meta_options = meta_df.columns.tolist()
        meta_options = [i for i in meta_options if i[1] in ['ORD','BIN','NOM']]
        meta_select = st.selectbox('meta', meta_options, key='tsne_meta')

    with st.expander("Data"):
        tabs = st.tabs(["Data", "Meta"])
        with tabs[0]:
            st.dataframe(data_df)
        with tabs[1]:
            st.dataframe(meta_df)
    
    with st.expander("Figure"):
        tsne_view(data_df, meta_df, meta_select, job_dir)





def app():
    data_dir = st.session_state.data_dir
    fasta_path = st.session_state.fasta_path
    # st.write(fasta_path)
    out_dir = st.session_state.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    job_dir = os.path.join(out_dir, 'dim_reduction')
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)
        

    gene_map = perform_get_gene_map(fasta_path)

    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)

    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
    else:
        st.title("Dimensionality Reduction")
        
        # protein_gene_map_path = readme[(readme['Class']=='Other_Map')&(readme['Data.Format']=="protein_gene_map")].iloc[0]['Path']
        protein_gene_map_path = next(
            iter(readme.loc[
                (readme['Class'] == 'Other_Map') & (readme['Data.Format'] == "protein_gene_map"), 
                'Path'
            ]), 
            "protein_gene_map.tsv"
        )
        
        protein_gene_map_path = Path(data_dir) / protein_gene_map_path
        
        if protein_gene_map_path.exists() and st.session_state.use_customized_gene_mapping:
            gene_map_table = pd.read_csv(protein_gene_map_path, sep='\t', header=None, names=['Protein','Gene'])
            gene_map = dict(zip(gene_map_table['Protein'], gene_map_table['Gene']))
            st.write(f"Using customized gene mapping from {protein_gene_map_path}")
        
            # st.dataframe(gene_map)
        
        tabs = st.tabs(["PCA",'T-SNE','UMAP'])
        with tabs[0]:
            single_pca(gene_map, readme, data_dir, job_dir)
            
        with tabs[1]:
            single_tsne(gene_map, readme, data_dir, job_dir)

        with tabs[2]:
            single_umap(gene_map, readme, data_dir, job_dir)