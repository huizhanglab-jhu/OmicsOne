import streamlit as st
import pandas as pd
import os,re,sys
import matplotlib.pyplot as plt
import seaborn as sns
from utils.process_readme import find_tn_pairs, union_genes_from_paths
import itertools as it
import numpy as np
from utils.omicsone_diff import compare_two_files
from utils.fasta import get_gene_map

def plot_heatmap_single(gene: str, diff_map: dict):
    data = []
    gene, gene_symbol = gene.split('@')
    for omics_select in sorted(diff_map.keys()):
        diff_df = diff_map[omics_select]
        diff_value = diff_df.loc[gene,:]['Log2FC(median)']
        row = [omics_select, diff_value]
        data.append(row)

    data = pd.DataFrame(data, columns=['Omics', 'Diff'])
    df = data.set_index('Omics')
    df = df.loc[:,['Diff']]
    # Create the heatmap within a specific figure
    fig, ax = plt.subplots(figsize=(1, 2))  # Set figure size
    sns.heatmap(df, annot=True, fmt=".2f", cmap="coolwarm", 
                cbar=False, ax=ax,  linewidths=0.5,  # Set line width for edges
                vmin=-1, vmax=1,
                linecolor='black')  # Set edge color))
    # Remove xlabel and xticklabels
    ax.set_xlabel("")  # Remove x-axis label
    ax.set_xticks([])  # Remove x-axis tick labels

    # Remove xlabel and xticklabels
    ax.set_ylabel("")  # Remove x-axis label
    ax.set_yticks([])  # Remove x-axis tick labels
    ax.set_title(gene_symbol)
    

    st.dataframe(data)
    # fig = sns.heatmap(data, x='Omics', y='Gene', z='Diff', cmap='coolwarm')
    return fig

def plot_heatmap_grid(genes: list, diff_map: dict, cols: int = 5):
    N = len(genes)
    rows = -(-N // cols)  # Calculate number of rows needed, round up
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 4))  # Grid size

    temp_tabs = st.tabs(genes)
    
    # Flatten axes for easier indexing
    axes = axes.flatten()
    
    for i in range(N):
        ax = axes[i]
        if i >= len(genes):
            axes[i].axis("off")  # Hide unused subplots
            continue
        
        gene_symbol = genes[i]

        data = []

        for omics_select in sorted(diff_map.keys()):
            diff_df = diff_map[omics_select]
            try:
                diff_value = diff_df[diff_df['Gene.Symbol'] == gene_symbol].iloc[0]['Log2FC(median)']
            except:
                diff_value = np.nan
            row = [omics_select, diff_value]
            data.append(row)
        
        data = pd.DataFrame(data, columns=['Omics', 'Diff'])
        df = data.set_index('Omics')
        df = df.loc[:,['Diff']]

        sns.heatmap(df, annot=True, fmt=".2f", cmap="coolwarm", 
                    cbar=False, ax=ax,  linewidths=0.5,  # Set line width for edges
                    vmin=-1, vmax=1,
                    linecolor='black')  
                # Remove xlabel and xticklabels
        ax.set_xlabel("")  # Remove x-axis label
        ax.set_xticks([])  # Remove x-axis tick labels

        # Remove xlabel and xticklabels
        ax.set_ylabel("")  # Remove x-axis label
        ax.set_yticks([])  # Remove x-axis tick labels
        ax.set_title(gene_symbol)

        with temp_tabs[i]:
            st.dataframe(data)
    
    plt.tight_layout()
    return fig


def heatmap_page(readme: pd.DataFrame,gene_map: dict,data_dir: str):
    # st.write(readme, gene_map, data_dir)
    tn_pairs = find_tn_pairs(readme)

    omics_options = sorted([i for i in tn_pairs.keys()])
    omics_selects = st.multiselect('omics', omics_options)

    if len(omics_selects) == 0:
        st.write('No omics selected!')
    else:
        output_dir = st.session_state['out_dir']
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        job_dir = os.path.join(output_dir,'omicsone_heatmap')
        if not os.path.exists(job_dir):
            os.makedirs(job_dir)

        # st.write(job_dir)

        gene_map = get_gene_map(st.session_state['fasta_path'])

        diff_map = dict()

        for omics_select in omics_selects:
            tumor_path = tn_pairs[omics_select]['Tumor']
            tumor_basename = os.path.splitext(tumor_path)[0]
            tumor_path = os.path.join(data_dir,tumor_path)
            normal_path = tn_pairs[omics_select]['Normal']
            normal_path = os.path.join(data_dir,normal_path)
            
            temp_diff_path = os.path.join(job_dir,f'{tumor_basename}_TNDiff.tsv')
            if os.path.exists(temp_diff_path):
                diff_df = pd.read_csv(temp_diff_path,sep='\t',index_col=0)
            else:
                diff_df = compare_two_files(tumor_path, normal_path, method = 'Wilcoxon(Unpaired)', fdr_cutoff=0.01, log2fc_cutoff=1)
                diff_df.to_csv(temp_diff_path,sep='\t')
            diff_map[omics_select] = diff_df

            if 'Gene.Symbol' not in diff_df.columns:
                diff_df['Gene.Symbol'] = [gene_map.get(i.split('.')[0])['gene'] if i.split('.')[0] in gene_map else 'NA' for i in diff_df.index]

        with st.expander("Data"):
            if len(omics_selects) > 0:
                tabs = st.tabs(omics_selects)
                for i,omics_select in enumerate(omics_selects):
                    with tabs[i]:
                        st.dataframe(diff_map[omics_select])

        paths = [[tn_pairs[i][t] for t in ["Tumor", "Normal"]] for i in omics_selects]
        paths =  list(it.chain.from_iterable(paths))
        paths = [os.path.join(data_dir,i) for i in paths if os.path.exists(os.path.join(data_dir,i))]
        # st.write(paths)

        genes = sorted(list(set(union_genes_from_paths(paths))))
        genes = [f"{i}@{gene_map.get(i.split('.')[0])['gene'] if i.split('.')[0] in gene_map else 'NA'}" for i in genes]
        # st.write(len(genes))

        tab_heatmap_single, tab_heatmap_multi = st.tabs(["Heatmap Single", "Heatmap Multi"])
        with tab_heatmap_single:
            gene = st.selectbox('gene', genes)
            fig = plot_heatmap_single(gene, diff_map)
            col_heatmap_single_fig, col_heatmap_single_settings = st.columns([1,5])
            with col_heatmap_single_fig:
                st.pyplot(fig)
            with col_heatmap_single_settings:
                pass
        with tab_heatmap_multi:
            gene_list_input = st.text_area("Enter your list of genes, separated by commas or tabs:", "")
            separator = st.selectbox('separator', [',',';', 'space'])

            if separator == 'space':
                genes = [str(i).strip() for i in re.split(r'\s+',gene_list_input)]
            else:
                genes = [str(i).strip() for i in gene_list_input.split(separator)]

            genes = [i for i in genes if i != ""]

            if len(genes) == 0:
                st.write("No genes selected!")
            else:
                fig = plot_heatmap_grid(genes, diff_map)
                st.pyplot(fig)


