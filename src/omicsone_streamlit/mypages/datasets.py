import streamlit as st
import os,sys, re
import pandas as pd
import numpy as np
# from utils.layout import display_user_projects
from config.update_settings import DATA_DIR, OUT_DIR, update_ini
from config.update_settings import FASTA_PATH, CHROM_PATH, CYTOBAND_PATH


def app():
    global DATA_DIR
    global OUT_DIR
    global FASTA_PATH
    global CYTOBAND_PATH
    global CHROM_PATH

    tab1, tab2, tab3, tab4 = st.tabs(["Input","Ouput","Reference","Misc"])


    with tab1:

        new_data_dir = st.text_input("Enter the new data directory path:")

        # Check if the directory exists and update settings
        if st.button("Update Dataset Folder"):
            if os.path.isdir(new_data_dir):
                DATA_DIR = new_data_dir
                update_ini("dirs","data_dir",new_data_dir)  # Persist the new path to settings.toml

                st.session_state.data_dir = new_data_dir

                st.success(f"Data directory updated to: {st.session_state.data_dir}")

            else:
                st.error(f"Invalid directory: {new_data_dir}")


        st.write(f"current data folder @ {DATA_DIR}")

        if os.path.isdir(DATA_DIR):
            readme_path = os.path.join(DATA_DIR,"readme.xlsx")
            if os.path.exists(readme_path):
                readme = pd.read_excel(readme_path)

                st.dataframe(readme, use_container_width=True)

    with tab2:

        out_dir = st.text_input("Enter the output directory path:")

        if st.button("Update Ouput Directory"):
            if os.path.isdir(out_dir):
                OUT_DIR = out_dir
                update_ini("dirs","out_dir",out_dir)  # Persist the new path to settings.toml
                st.session_state.out_dir = out_dir

                st.success(f"Data directory updated to: {st.session_state.out_dir}")

            else:
                st.error(f"Invalid directory: {out_dir}")


        st.write(f"current data folder @ {OUT_DIR}")
    
    with tab3:

        fasta_path = st.text_input("Enter the protein database path: (.fasta)")

        if st.button("Update Protein Path"):
            if os.path.isfile(fasta_path):
                FASTA_PATH = fasta_path
                update_ini("paths","fasta_path",fasta_path)  # Persist the new path to settings.toml
                st.session_state.fasta_path = fasta_path

                st.success(f"Protein database updated to: {st.session_state.fasta_path}")

            else:
                st.error(f"Invalid path: {fasta_path}")
        st.write(f"current fasta_path @ {FASTA_PATH}")

        chrom_path = st.text_input("Enter the chromosomes path: (chromosomes.txt)")

        if st.button("Update Chromosomes Path"):
            if os.path.isfile(chrom_path):
                CHROM_PATH = chrom_path
                update_ini('paths',"chrom_path",chrom_path)
                st.session_state.chrom_path = chrom_path
                st.success(f"The chromosomes path updated to: {st.session_state.chrom_path}")
            else:
                st.error(f"Invalid path: {chrom_path}")
        st.write(f"current chromosomes path @ {CHROM_PATH}")

        cytoband_path = st.text_input("Enter the cytoband path: (cytoband.txt)")

        if st.button("Update Cytoband Path"):
            if os.path.isfile(cytoband_path):
                CYTOBAND_PATH = cytoband_path
                update_ini('paths',"cytoband_path",cytoband_path)
                st.session_state.cytoband_path = cytoband_path
                st.success(f"The cytoband path updated to: {st.session_state.cytoband_path}")
            else:
                st.erro(f"Invalid path: {cytoband_path}")
        st.write(f"current cytoband path @ {CYTOBAND_PATH}")
        
    with tab4:
        use_customized_gene_mapping = st.checkbox("Use customized protein_gene_map.tsv if available in the folder", value=True)
        update_ini('misc',"use_customized_gene_mapping", str(use_customized_gene_mapping))
        st.session_state.use_customized_gene_mapping = use_customized_gene_mapping
        








  