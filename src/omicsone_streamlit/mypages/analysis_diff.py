import streamlit as st
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils.fasta import get_gene_map
from utils.omicsone_boxplot import boxplot_page
from utils.omicsone_heatmap import heatmap_page
from utils.omicsone_volcano import volcano_page


@st.cache_data
def perform_get_gene_map(path):
    gene_map = get_gene_map(path)
    return gene_map


def app():


    data_dir = st.session_state.data_dir
    fasta_path = st.session_state.fasta_path
    # st.write(fasta_path)

    gene_map = perform_get_gene_map(fasta_path)
    
    # st.write(gene_map)

    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)

    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
        st.stop()

    if st.session_state.use_customized_gene_mapping:
        # protein_gene_map_path = readme[(readme['Class']=='Other_Map')&(readme['Data.Format']=="protein_gene_map")].iloc[0]['Path']
        protein_gene_map_path = next(
            iter(readme.loc[
                (readme['Class'] == 'Other_Map') & (readme['Data.Format'] == "protein_gene_map"), 
                'Path'
            ]), 
            "protein_gene_map.tsv"
        )
        protein_gene_map_path = os.path.join(data_dir, protein_gene_map_path)
        if os.path.exists(protein_gene_map_path):
            gene_map_table = pd.read_csv(protein_gene_map_path, sep='\t', header=None, names=['Protein','Gene'])
            gene_map = dict()
            for protein, gene in zip(gene_map_table['Protein'], gene_map_table['Gene']):
                gene_map[protein] = {
                    'gene': gene
                }
            # gene_map = dict(zip(gene_map_table['Protein'], gene_map_table['Gene']))
            st.write(f"Using customized gene mapping from {protein_gene_map_path}")
            
    tab_boxplot, tab_heatmap, tab_volcano = st.tabs(["Boxplot", "Heatmap","Volcano"])
    with tab_boxplot:
        boxplot_page(readme,gene_map,data_dir)
    with tab_heatmap:
        heatmap_page(readme,gene_map,data_dir)
    with tab_volcano:   
        volcano_page(readme, gene_map, data_dir)
        


def mixed_app():
    data_dir = st.session_state.data_dir
    fasta_path = st.session_state.fasta_path
    # st.write(fasta_path)

    gene_map = perform_get_gene_map(fasta_path)

    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)

    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
        st.stop()

    tab_boxplot, tab_heatmap, tab_volcano = st.tabs(["Boxplot", "Heatmap","Volcano"])
    # with tab_boxplot:
    #     boxplot_page(readme,gene_map,data_dir)
    # with tab_heatmap:
    #     heatmap_page(readme,gene_map,data_dir)
    with tab_volcano:   
        volcano_page(readme, gene_map, data_dir)