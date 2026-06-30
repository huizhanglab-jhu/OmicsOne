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
from utils.timestamp import build_timestamp
import configparser
import subprocess

@st.cache_data
def perform_compare_two_files(file_a,file_b, params):
    return compare_two_files(file_a, file_b, method = params['method'], 
                      fdr_cutoff=params['fdr_cutoff'], log2fc_cutoff=params['log2fc_cutoff'])
    

def volcano_plot(diff: pd.DataFrame, params : dict):
    height = params.get('height',4)
    width = params.get('width',4)
    dpi = params.get('dpi',200)
    colors = params.get('colors',['grey','red','blue'])
    fig, ax = plt.subplots(figsize=(width, height),dpi=dpi)
    sns.set_style('white')
    up_count = diff[diff['Significance']=='S-U'].shape[0]
    down_count = diff[diff['Significance']=='S-D'].shape[0]
    color0 = colors[0]
    sns.scatterplot(data=diff, x='Log2FC(median)', y='-Log10(FDR)',
                    color=color0, s=1)
    color1 = colors[1]
    color2 = colors[2]
    sns.scatterplot(data=diff[diff['Significance']=='S-U'], x= 'Log2FC(median)', y='-Log10(FDR)',
                    color=color1,s=5,label='UP({})'.format(up_count))
    sns.scatterplot(data=diff[diff['Significance']=='S-D'], x= 'Log2FC(median)', y='-Log10(FDR)',
                color=color2,s=5,label='DOWN({})'.format(down_count))
    
    x_upper = 1.1 * np.max([abs(np.nanmin(diff['Log2FC(median)'])), np.nanmax(diff['Log2FC(median)'])])
    x_down = -1 * x_upper
    y_upper = np.nanmax(diff['-Log10(FDR)']) * 1.1
    
    plt.plot([0, 0], [0, y_upper], color='k', linewidth=0.5, linestyle="--")
    plt.plot([-1, -1], [0, y_upper], color='k', linewidth=0.5, linestyle="--")
    plt.plot([1, 1], [0, y_upper], color='k', linewidth=0.5, linestyle="--")
    # plt.plot([-5, 5], [0, 0], color='k', linewidth=0.5, linestyle="--")
    plt.plot([-5, 5], [-np.log10(0.01), -np.log10(0.01)], color='k', linewidth=0.5, linestyle="--")
    plt.legend()
    title = params.get('title','differential expression analysis')
    plt.title(title)
    xlabel = params.get('xlabel', 'Log2FC(Tumor/NAT)')
    plt.xlabel(xlabel)

    return fig

@st.cache_data
def perform_enrichment(pure_up_genes, pure_down_genes, total_genes, job_dir):
    import time
    from utils.pathway import omicsone_enrichr
    up_enrichr_df = omicsone_enrichr(pure_up_genes,total_genes,job_dir,gene_sets=["MSigDB_Hallmark_2020"])
    time.sleep(5)
    down_enrichr_df = omicsone_enrichr(pure_down_genes,total_genes,job_dir,gene_sets=['MSigDB_Hallmark_2020'])
    return up_enrichr_df, down_enrichr_df

@st.cache_data
def perform_get_gene_map(path):
    gene_map = get_gene_map(path)
    return gene_map

def volcano_page(readme: pd.DataFrame,gene_map: dict,data_dir: str):

    out_dir = st.session_state.out_dir
    job_dir = os.path.join(out_dir, "diff_volcano")
    
    fasta_path = st.session_state.fasta_path
    
    # gene_map = perform_get_gene_map(fasta_path)
    
    # if st.session_state.use_customized_gene_mapping:
    #     protein_gene_map_path = readme[(readme['Class']=='Other_Map')&(readme['Data.Format']=="protein_gene_map")].iloc[0]['Path']
    #     protein_gene_map_path = os.path.join(data_dir, protein_gene_map_path)
    #     if os.path.exists(protein_gene_map_path):
    #         gene_map_table = pd.read_csv(protein_gene_map_path, sep='\t', header=None, names=['Protein','Gene'])
    #         gene_map = dict(zip(gene_map_table['Protein'], gene_map_table['Gene']))
    #         st.write(f"Using customized gene mapping from {protein_gene_map_path}")
    
    # st.write(gene_map)

    if "volcano_colors" not in st.session_state:
        st.session_state.volcano_colors = ["#808080","#FF0000", "#0000FF"]

    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    omics_options = [i for i in set(readme['Class']) if not re.search('Other.',i)]
    omics_select = st.selectbox('omics',omics_options, key='diff.volcano')    

    col_a, col_b = st.columns(2)
    
    group_a_options = list(readme[readme['Class']==omics_select]['Path'])

    with col_a:
        group_a_select = st.selectbox('A',group_a_options, key='diff.volcano.a') 

    row_a = readme[readme['Path']==group_a_select].iloc[0]
    a_Class = row_a['Class']
    a_EM = row_a['Experiment.Method']
    a_Group = row_a['Group']	
    a_QM = row_a['Quant.Method']
    a_log = row_a['logTransform']
    a_Pathology =  row_a['Pathology']

    rows_b = list(readme[(readme['Class']==a_Class)&(readme['Experiment.Method']==a_EM)&(
        readme['Group'] == a_Group)&(readme['Quant.Method']==a_QM) & (readme['logTransform']==a_log
    )]['Path'])
    
    group_b_options = [i for i in rows_b if i != group_a_select]
    with col_b:
        group_b_select = st.selectbox('B',group_b_options, key='diff.volcano.b')
    row_b = readme[readme['Path']==group_b_select].iloc[0]
    b_Pathology = row_b["Pathology"]

    file_a  = os.path.join(data_dir, group_a_select)

    file_b = os.path.join(data_dir, group_b_select)

    method_options = ['Wilcoxon(Unpaired)','Wilcoxon(Paired)','T-test(Unpaired)','T-test(Paired)',]
    with st.expander("Settings"):
        method  = st.selectbox(label="Method", options  = method_options)
        fdr_cutoff = st.number_input(label="FDR threshold", value = 0.01, step=0.01 )
        log2fc_cutoff = st.number_input(label = "Log2 fold change threshold", value = 1.0, step=0.1)

    timestamp = build_timestamp()

    diff_params = {
        "method" : method,
        "fdr_cutoff" : fdr_cutoff,
        "log2fc_cutoff": log2fc_cutoff

    }

    diff_path = os.path.join(job_dir, f"diff_{timestamp}.tsv")    
    diff_log_path = os.path.join(job_dir, f"diff_log_{timestamp}.ini")

    # Create a ConfigParser object
    config = configparser.ConfigParser()
    config["Parameters"] = diff_params
    # Write parameters to the .ini file
    with open(diff_log_path, "w") as configfile:
        config.write(configfile)


    with st.spinner("Processing... Please wait."):
        
        data_a = pd.read_csv(file_a, sep="\t", header=0, index_col=0)
        data_a.columns = [i + '.A' for i in data_a.columns.values]
        data_b = pd.read_csv(file_b, sep="\t", header=0, index_col=0)
        data_b.columns = [i + '.B' for i in data_b.columns.values]
        data = pd.concat([data_a, data_b], axis=1)
        
    with st.expander("Data"):
        tabs = st.tabs(["A","B","Combined"])
        with tabs[0]:
            st.dataframe(data_a)
        with tabs[1]:
            st.dataframe(data_b)
        with tabs[2]:
            st.dataframe(data)
            
    diff_df = perform_compare_two_files(file_a, file_b, diff_params)
    diff_df.to_csv(diff_path,sep='\t')
    st.write(f"DE analysis has been saved to {diff_path}")

    with st.expander("Results"):
        st.dataframe(diff_df)

    with st.expander("Figure"):
        cols = st.columns([1,1])
        with cols[1]:
            title = st.text_input("title", "differential expression analysis ")
            figure_width = st.number_input(label="figure_width", value=4, step=1)
            figure_height = st.number_input(label = "figure_height",value = 4, step=1)
            st.session_state.volcano_colors[0] = st.color_picker(label="color.background", value=st.session_state.volcano_colors[0])
            st.session_state.volcano_colors[1] = st.color_picker(label="color.up", value=st.session_state.volcano_colors[1])
            st.session_state.volcano_colors[2] = st.color_picker(label="color.down", value=st.session_state.volcano_colors[2])
            if st.button('reset colors',key='volcano_default_colors'):
                st.session_state.volcano_colors = ["#808080","#FF0000", "#0000FF"]


        with cols[0]:
            params = {
                'title': title,
                "width": figure_width,
                "height": figure_height,
                "colors": st.session_state.volcano_colors
            }
            fig = volcano_plot(diff_df,params)
            st.pyplot(fig)

    with st.expander("Enrichment"):
        tabs = st.tabs(['GSEApy'])
        with tabs[0]:
            st.write('gseapy')
            
            diff = diff_df.copy()
            total_genes = set([i for i in data.index])
            up_genes = set([i for i in diff[diff['Significance']=='S-U'].index])
            down_genes = set([i for i in diff[diff['Significance']=='S-D'].index])
            pure_down_genes = down_genes - up_genes
            pure_up_genes = up_genes - down_genes
            
            pure_down_genes = [gene_map[i.split(".")[0]]["gene"] for i in pure_down_genes if i.split(".")[0] in gene_map]
            pure_up_genes = [gene_map[i.split(".")[0]]["gene"] for i in pure_up_genes if i.split(".")[0] in gene_map]
            total_genes = [gene_map[i.split(".")[0]]["gene"] for i in total_genes if i.split(".")[0] in gene_map]
            # st.write("Total genes: ", list(total_genes)[:4])
            # st.write(len(total_genes),len(up_genes),len(down_genes),len(pure_down_genes),len(pure_up_genes))
            # st.write(pure_down_genes[:4])
            # st.write(pure_up_genes[:4])
            up_enrichr_df, down_enrichr_df = perform_enrichment(pure_up_genes, pure_down_genes, total_genes, job_dir)
            tabs = st.tabs(['Up','Down',"Plot"])
            with tabs[0]:
                st.dataframe(up_enrichr_df)
            with tabs[1]:
                st.dataframe(down_enrichr_df)
            with tabs[2]:
                from utils.pathway import plot_enrichr_both
                skip_pathways = [
                    'Phagosome','Human papillomavirus infection','Pertussis',
                    'Malaria','Arrhythmogenic right ventricular cardiomyopathy',
                    'Staphylococcus aureus infection','Regulation of actin cytoskeleton'
                ]
                up_enrichr_df2 = up_enrichr_df[~up_enrichr_df['Term'].isin(skip_pathways)].sort_values('Adjusted P-value').head(10)
                down_enrichr_df2 = down_enrichr_df[~down_enrichr_df['Term'].isin(skip_pathways)].sort_values('Adjusted P-value').head(10)
                fig = plot_enrichr_both(up_enrichr_df2, down_enrichr_df2, title="Enrichment Analysis")

                st.pyplot(fig)

        # with tabs[0]:
        #     project_name = st.text_input(label="Project Name", value="WebGestAltR")
        #     gene_set_options = ["Sig-Up","Up","Sig-Down","Down"]
        #     gene_set_option = st.selectbox(label="Gene Set", options = gene_set_options)

        #     function_tabs = st.tabs(['Over-Representative Analysis',
        #                              'Gene Set Enrichment Analysis',
        #                              'Network Topology-based Analysis'])
            


        #     with function_tabs[0]:
        #         working_dir = os.path.join(job_dir, "webgestatltr_ora")
        #         if not os.path.exists(working_dir):
        #             os.mkdir(working_dir)

        #         gene_list = []
        #         reference_list = []
        #         gene_list_path  = os.path.join(working_dir,"interesting_gene_list.txt")
        #         reference_list_path = os.path.join(working_dir, "reference_gene_list.txt")
        #         if gene_set_option == "Sig-Up":
        #             gene_list = list(diff_df[diff_df['Significance']=='S-U'].index)
        #         elif gene_set_option == "Sig-Down":
        #             gene_list = list(diff_df[diff_df["Significance"]=="S-D"].index)
        #         elif gene_set_option == "Up":
        #             gene_list = list(diff_df[(diff_df["Significance"]=="U")|(diff_df["Significance"]=="S-U")].index)
        #         elif gene_set_option == "Down":
        #             gene_list = list(diff_df[(diff_df["Significance"]=="D")|(diff_df["Significance"]=="S-D")].index)
        #         with open(gene_list_path,"w") as f:
        #             for gene in gene_list:
        #                 f.write(f'{gene.split(".")[0]}\n')
        #             f.close()

        #         reference_list = list(diff_df.index)
        #         with open(reference_list_path,"w") as f:
        #             for gene in reference_list:
        #                 f.write(f'{gene.split(".")[0]}\n')
        #             f.close()

        #         # Define the command and its arguments
        #         r_script_dir = r'/Users/yingweihu/project/webgestaltR'
        #         r_script_path = os.path.join(r_script_dir,'run_webgestalt.R')
                
        #         rscript_command = [
        #             "Rscript",  # Command to run R script
        #             r_script_path,  # Path to your R script
        #             gene_list_path,  # Path to gene file
        #             reference_list_path,  # Path to reference file
        #             working_dir,  # Output directory
        #             project_name # Project name
        #         ]

        #         ora_tabs = st.tabs(['Log','Report','GO'])

        #         with ora_tabs[0]:
                        
        #             # Run the command
        #             try:
        #                 result = subprocess.run(
        #                     rscript_command,
        #                     capture_output=True,  # Capture stdout and stderr
        #                     text=True  # Get output as string
        #                 )

        #                 # Print the output
        #                 st.write("Standard Output:\n", result.stdout)
        #                 st.write("Standard Error:\n", result.stderr)
                    
        #                 # Check if the command was successful
        #                 if result.returncode == 0:
        #                     st.write("R script executed successfully.")
        #                 else:
        #                     st.write("R script execution failed.")
        #             except FileNotFoundError:
        #                 print("Error: Rscript not found. Ensure that Rscript is installed and available in your PATH.")


        #         with ora_tabs[1]:
                    
        #             html_file_path = os.path.join(working_dir, f"Project_{project_name}", f"Report_{project_name}.html")
        #             if os.path.exists(html_file_path):
        #                 # Read the content of the HTML file
        #                 with open(html_file_path, "r", encoding="utf-8") as html_file:
        #                     html_content = html_file.read()

        #                 # Display the HTML content in Streamlit
        #                 st.components.v1.html(html_content, height=1000,scrolling=True)  # Adjust the


        #         with ora_tabs[2]:
                    
        #             go_summary_path = os.path.join(working_dir, f"Project_{project_name}", f"goslim_summary_{project_name}.png")
        #             st.image(go_summary_path)

        #     with function_tabs[1]:
        #         pass

        #     with function_tabs[2]:
        #         pass

    st.write("END")

    
