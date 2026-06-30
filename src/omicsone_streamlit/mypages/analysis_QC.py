import streamlit as st
from utils.fasta import get_gene_map
import os, re,sys
import pandas as pd
import numpy as np
from utils.fasta import get_gene_map
from myviews.qc_views import qc_corr_page
from myviews.qc_views import statistics_view as qc_stats
from myviews.qc_views import missing_values_view as qc_missing_values
from myviews.qc_views import reproducibility_view as qc_reproducibility
from myviews.qc_views import pca_view as qc_pca
from myviews.qc_views import QCInput
from pathlib import Path

@st.cache_data
def perform_get_gene_map(fasta_path):
    gene_map = get_gene_map(fasta_path)
    return gene_map

@st.cache_data
def read_readme(data_dir):
    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)
    else:
        print(f"Data directory {data_dir} does not exist.")
        
    return readme


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

    # gene_map = perform_get_gene_map(fasta_path)

    # readme = None
    # if os.path.isdir(data_dir):
    #     readme_path = os.path.join(data_dir,"readme.xlsx")
    #     if os.path.exists(readme_path):
    #         readme = pd.read_excel(readme_path)

    readme = read_readme(data_dir)
    
    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
        st.stop()
        
        
    st.title("Quaity Control")
    
    file_cols = st.columns([1,2])
    with file_cols[0]:
        omics_list = [i for i in readme['Class'].unique() if not re.search('^Other',i)]
        omics_select = st.selectbox('Omics.Type', omics_list)

    with file_cols[1]:
        file_list = readme[readme['Class']==omics_select]['Path'].unique()
        file_select = st.selectbox('Data', file_list)
        
    meta_cols = st.columns([1,1,1])
    with meta_cols[0]:
        meta_list = readme[readme['Class']=='Other_Meta']['Path'].unique()
        sel_meta = st.selectbox('Meta', meta_list)
        
    meta_full_path = Path(data_dir) / sel_meta
    meta = pd.read_csv(str(meta_full_path), sep='\t',index_col=0, header=[0,1])
    
    set_options = [f"{i}@{j}" for i,j in  meta.columns if j in ["ORD","BIN","NOM"]]
    with meta_cols[1]:
        sel_set = st.selectbox('Set', set_options)
        
    group_options = [f"{i}@{j}" for i,j in  meta.columns if j in ["ORD","BIN","NOM"]]
    with meta_cols[2]:
        sel_group = st.selectbox('Group', group_options)


    sel_group = tuple(sel_group.split('@'))
    sel_set = tuple(sel_set.split('@'))

    path = os.path.join(data_dir, file_select)
    
    pathology_type = readme[readme['Path']==file_select]['Pathology'].values[0]

    # data = pd.read_csv(path,sep="\t",index_col=0)
    
    if omics_select in ["Glycoproteomics"]:
        data = pd.read_csv(path,sep="\t",index_col=[0,1])
    elif omics_select in ["Phosphoproteomics"]:
        data = pd.read_csv(path,sep="\t",index_col=[0,1,2,3])
    else:
        data = pd.read_csv(path,sep="\t",index_col=0)
    
    data.columns = [str(i).strip() for i in data.columns]

    # data2 = data.iloc[:,:3]
    samples = list(data.columns.values)

    with st.expander("Data"):
        data_tabs = st.tabs(["Data","Meta"])
        # from utils.display import display_table
        # from st_aggrid import AgGrid
        with data_tabs[0]:
            st.dataframe(data)
            pass
            # AgGrid(data)
        with data_tabs[1]:
            st.dataframe(meta)
            
    options = ["SampleCorrelation","Statistics","MissingValues","Reproducibility","PCA"]
    cols = st.columns(len(options))
    
    
    tab_names = []
    for i, option in enumerate(options):
        with cols[i]:
            show_stats = st.checkbox(option, value = st.session_state[f'QC_{option}'])
            if show_stats != st.session_state[f'QC_{option}']:
                st.session_state[f'QC_{option}'] = show_stats
            if st.session_state[f'QC_{option}']:
                tab_names.append(option)
                
    tabs = []
    if len(tab_names) > 0:
        tabs = st.tabs(tab_names)
    
    cfg = QCInput(
                omics_select=omics_select,
                data=data,
                meta=meta,
                sel_set=sel_set,
                sel_group=sel_group,
                pathology_type=pathology_type
            )
                
    if "SampleCorrelation" in tab_names:
        with tabs[tab_names.index("SampleCorrelation")]:
            qc_corr_page( samples,  data)
    
    if "Statistics" in tab_names:
        with tabs[tab_names.index("Statistics")]:
            qc_stats(cfg)
            
    if "MissingValues" in tab_names:
        with tabs[tab_names.index("MissingValues")]:
            qc_missing_values(cfg)
        
    if "Reproducibility" in tab_names:
        with tabs[tab_names.index("Reproducibility")]:
            qc_reproducibility(cfg)
        
    if "PCA" in tab_names:
        with tabs[tab_names.index("PCA")]:
            qc_pca(cfg)
