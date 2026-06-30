import streamlit as st
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils.omicsx_gene_corr import calculate_corr, plot_corr, gsea_prerank
from utils.fasta import get_gene_map
from utils.diff import compare_two_groups
from utils.omicsx_group_diff import get_matching_rows, get_matching_path_with_pathology
from adjustText import adjust_text
from utils.omicsx_commons import read_omics2
from utils.omicsx_sample_corr import calculate_sample_corr, calculate_gene_sample_corr
import numpy as np
from matplotlib.colors import ListedColormap

# import gseapy
@st.cache_data
def perform_two_diff_analysis(omics1_merged, omics2_merged, samples_a_a, samples_a_b, 
                         samples_b_a, samples_b_b):
    """
    Perform differential analysis and create visualization
    """
    # calculate statistic test
    omics1_diff = compare_two_groups(omics1_merged, samples_a_a, samples_a_b, 
                                   method="Wilcoxon(Unpaired)",
                                   max_miss_ratio_global=0.5, 
                                   max_miss_ratio_group=0.5,
                                   fdr_cutoff=0.01, 
                                   log2fc_cutoff=1)

    omics2_diff = compare_two_groups(omics2_merged, samples_b_a, samples_b_b,
                                   method="Wilcoxon(Unpaired)",
                                   max_miss_ratio_global=0.5, 
                                   max_miss_ratio_group=0.5,
                                   fdr_cutoff=0.01, 
                                   log2fc_cutoff=1)

    return omics1_diff, omics2_diff

@st.cache_data
def perform_corr_prank(path1, path2, gene_map, out_dir):
    df = calculate_corr(path1, path2, out_dir=out_dir, gene_map=gene_map)
    job_dir = os.path.join(out_dir,"genecorrelation","gsea1")
    gsea_prerank_df = gsea_prerank(df,col_name="Gene Correlation",out_dir=job_dir)
    return df, gsea_prerank_df



@st.cache_data
def plot_group_diff(merged_diff, class_a_select, class_b_select, group_diff_select1, group_diff_select2):
                   # Create scatter plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot log2FC values
    ax.scatter(merged_diff['Log2FC(median)_omics1'], 
            merged_diff['Log2FC(median)_omics2'],
            alpha=0.5)

    # Add labels and title
    ax.set_xlabel(f'log2FC {class_a_select}')
    ax.set_ylabel(f'log2FC {class_b_select}')
    ax.set_title(f'Differential Analysis: {group_diff_select1} vs {group_diff_select2}')

    # Add a horizontal and vertical line at y=0 and x=0
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)

    # Optional: Add feature labels for significant points
    # You can adjust the thresholds as needed
    significant_mask = (abs(merged_diff['Log2FoldChange']) > 5) 
    
    texts = []
    
    for idx, row in merged_diff[significant_mask].iterrows():
        x = row['Log2FC(median)_omics1']
        y = row['Log2FC(median)_omics2']
        texts.append(ax.text(x, y, row['Gene'], fontsize=8))

    # Adjust text positions to minimize overlaps
    adjust_text(texts, 
               arrowprops=dict(arrowstyle='->', color='gray', lw=0.5),
               expand_points=(1.5, 1.5),
               force_points=(0.5, 0.5))

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    return fig

def perform_get_gene_map(fasta_path):
    gene_map = get_gene_map(fasta_path=fasta_path)
    return gene_map


def perform_add_gene_name(df, gene_map, column_name=None):
    if column_name is None:
        df['Gene'] = [gene_map.get(i.split(".")[0],
                    {'gene':'NA'})['gene'] for i in df.index]
    elif column_name in df.columns:
        df[column_name] = [gene_map.get(i.split(".")[0],
                    {'gene':'NA'})['gene'] for i in df[column_name]]
    else:
        st.write(f"Column {column_name} not found in dataframe, Failed to add gene name")
    return df

@st.cache_data
def perform_sample_corr(omics1, omics2, out_dir):
    df_samplecorr = calculate_sample_corr(omics1, omics2, out_dir)
    return df_samplecorr

@st.cache_data
def perform_calculate_gene_sample_corr(sample_corr_df, omics1, omics2, out_dir):
    return calculate_gene_sample_corr(sample_corr_df, omics1, omics2, out_dir)


def plot_sample_corr(df_samplecorr, ann_df, cols=['Sex','Stage']):
        # Create a figure with a specific size
    fig = plt.figure(figsize=(10, 6))
    n_cols = len(cols)
    # Create two axes with specific height ratios
    # The first axis (ax1) will be larger for the bar plot
    # The second axis (ax2) will be smaller for the heatmap
    gs = fig.add_gridspec(n_cols+1, 1, 
        height_ratios=[20] + [1]*n_cols, 
    hspace=0.02)
    ax1 = fig.add_subplot(gs[0])


    # Example data (replace with your actual data)
    # For bar plot
    bar_data = df_samplecorr.sort_values('Corr', ascending=False)

    # For heatmap
    # heatmap_data = pd.DataFrame(
    #     np.random.choice(['Status1', 'Status2', 'Status3'], size=(2, 50)),
    #     index=['Survival', 'PlatinumStatus'],
    #     columns=bar_data.index
    # )

    # Create bar plot
    sns.barplot(data=bar_data, x=bar_data.index, y='Corr', 
                color='blue', ax=ax1)

    # Remove x-axis labels from bar plot
    ax1.set_xticklabels([])
    ax1.set_xlabel('')

    # Customize bar plot
    ax1.set_title('Sample-wise mRNA-protein correlation')
    ax1.set_ylabel('Index (Spearman correlation)')

    legends_info = []
    for i, col in enumerate(cols):
        ax2 = fig.add_subplot(gs[i+1])
        # Select and reorder annotation data to match bar plot order
        ann_df2 = ann_df[cols].loc[bar_data.index]
 
        unique_values = ann_df2[col].unique()
        # Create a categorical color palette
        palette = sns.color_palette("Set3", n_colors=len(unique_values))
        # Create mapping dictionary
        color_dict = dict(zip(unique_values, palette))
        
        # Convert categorical values to their corresponding colors
        color_matrix = np.array([color_dict[val] for val in ann_df2[col]]).reshape(1, -1)
        
        # Plot this row of the heatmap
        sns.heatmap(color_matrix, 
                    cmap=ListedColormap(palette),
                    cbar=False,
                    # cbar_kws={'label': col},
                    yticklabels=[col],
                    xticklabels=False,
                    ax=ax2)
        ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0)
        
        # Create custom legend
        from matplotlib.patches import Patch
        # legend_elements = [Patch(facecolor=color_dict[val], 
        #                         label=val) for val in unique_values]
        # Store legend info for this annotation
        legend_elements = [Patch(facecolor=color_dict[val], label=val) for val in unique_values]
        legends_info.append((col, legend_elements))

    legend_y = 0.25
    for col, legend_elements in legends_info:
        leg =fig.legend(handles=legend_elements, 
                   title=col,
                   loc='lower left',
                   alignment='left',
                   bbox_to_anchor=(0.15, legend_y),
                   ncol=len(legend_elements),
                   frameon=False)

        legend_y -= 0.07  # Adjust spacing between legend lines as needed

        
    plt.subplots_adjust(bottom=0.35)  # Make room for legends at the top
        
    # Adjust layout
    plt.tight_layout()
    return fig

@st.cache_data
def perform_merge_col_colors(my_col_colors1, my_col_colors2, my_col_colors3, my_col_colors4):
    # Add suffixes to distinguish columns from different merges
    types = pd.merge(my_col_colors1, my_col_colors2, how='left', left_index=True, right_index=True, 
                    suffixes=('_1', '_2'))
    types = pd.merge(types, my_col_colors3, how='left', left_index=True, right_index=True,
                    suffixes=('', '_3'))
    types = pd.merge(types, my_col_colors4, how='left', left_index=True, right_index=True,
                    suffixes=('', '_4'))
    types.columns=["group1","group2","group3","group4"]
    # if len(sys.argv)>3:
    #     types = pd.merge(types, annotation.T, how='left', left_index=True, right_index=True)
    return types

@st.cache_data
def plot_merge_col_colors(types, out_dir):
    import matplotlib
        # First create a mapping of unique values to numbers
    unique_values = np.unique(types.values)
    value_to_num = {val: i for i, val in enumerate(unique_values)}

    # Convert the data to numerical values
    types_numeric = types.replace(value_to_num)

    # Create custom colormap from the original colors
    colors = ['green', 'red', 'cyan', 'magenta', 'yellow', 'blue', 'black']
    custom_cmap = matplotlib.colors.ListedColormap(colors[:len(unique_values)])

    # Create the plot
    fig, ax3 = plt.subplots(figsize=(5,2))
    plt.tight_layout()

    sns.heatmap(types_numeric.T, 
                cmap=custom_cmap,
                cbar=False,
                xticklabels=False,
                linewidths=0.1,
                linecolor='black')
    ax3.set_xlabel('')

    ax3.set_yticklabels(types_numeric.columns, rotation=0, fontsize=10)
    plt.savefig(os.path.join(out_dir,"figure7-sns.png"), dpi=100, bbox_inches='tight')
    return fig

@st.cache_data
def perform_plot_cluster_comparison(types, out_dir):
    from sklearn.metrics.cluster import adjusted_rand_score
    df_ari = pd.DataFrame(columns=types.columns,index=types.columns)
    for each in types.columns:
        for each2 in types.columns:
            df_ari.at[each,each2]=adjusted_rand_score(types[each],types[each2])
    df_ari=df_ari.astype(float).round(2)
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.heatmap(df_ari.astype(float),vmax=1,vmin=-1,cmap='RdBu_r', annot=True,
                cbar_kws={'label': 'Adjusted Rand Index'}
                    )
    plt.savefig(os.path.join(out_dir, "ARI.png"),dpi=100,bbox_inches='tight')
    return fig, df_ari

def app():

    data_dir = st.session_state.data_dir
    out_dir = st.session_state.out_dir

    fasta_path = st.session_state.fasta_path

    gene_map = perform_get_gene_map(fasta_path)
    

    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)
    

    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
    else:

        col_a, col_b = st.columns(2)


        with col_a:
            class_a_options = [i for i in set(readme['Class']) if not re.search('Other.',i)]
            class_a_select = st.selectbox('omics-1',class_a_options)  

            group_a_options = list(readme[readme['Class']==class_a_select]['Path'])
            path_a_select = st.selectbox('file-A',group_a_options)    

        with col_b:     
            class_b_options = [i for i in class_a_options if i != class_a_select]
            class_b_select = st.selectbox('omics-2',class_b_options)

            group_b_options = list(readme[readme['Class']==class_b_select]['Path'])
            path_b_select = st.selectbox('file-B', group_b_options)

        annotation_options = list(readme[readme['Class'].str.contains('Other')]['Path'])
        annotation_options = [i for i in annotation_options if not re.search('CaseList',i)]
        annotation_select = st.selectbox('Annotation',annotation_options)

        annotation_path = os.path.join(data_dir, annotation_select)
        ann_df = pd.read_csv(annotation_path,sep="\t", header=[0,1])

        # Extract column names and data types
        columns = ann_df.columns.get_level_values(0)  # First row
        data_types = ann_df.columns.get_level_values(1)  # Second row

        # Rename columns to use only the names
        ann_df.columns = columns

        # Store data types in a dictionary for later use
        data_type_dict = dict(zip(columns, data_types))

        ann_df = ann_df.set_index('case_id')

        path1 = os.path.join(data_dir,path_a_select)
        path2 = os.path.join(data_dir,path_b_select)

        omicsx_commons_dir = os.path.join(out_dir,"omicsx_commons")
        if not os.path.exists(omicsx_commons_dir):
            os.mkdir(omicsx_commons_dir)
        omics1, omics2 = read_omics2(path1, path2, out_dir=omicsx_commons_dir)


        # Create tabs
        tab1, tab2, tab3, tab4 = st.tabs(["Gene Correlation", "Sample Correlation", 
                                    "Sample Clustering", "Group Difference"])

        # Content for each tab
        with tab1:
            gene_corr_dir = os.path.join(out_dir,"omicsx_gene_correlation")
            if not os.path.exists(gene_corr_dir):
                os.mkdir(gene_corr_dir)
            

            df, gsea_prerank_df = perform_corr_prank(path1, path2, gene_map, gene_corr_dir)

            col1 , col2 = st.columns(2)

            with col1:
                fig = plot_corr(df, out_dir=gene_corr_dir)

                st.pyplot(fig)

            with col2:
                st.write("Gene-wise correlation")
                df["P"] = df["P"].apply(lambda x: f"{x:.2e}")
                df["BH adjusted P"] = df["BH adjusted P"].apply(lambda x: f"{x:.2e}")
                st.dataframe(df)

            st.write("Enriched pathways related to Gene-wise correlation")

            st.dataframe(gsea_prerank_df.res2d)



        with tab2:
            sample_corr_dir = os.path.join(out_dir,"omicsx_sample_correlation")
            if not os.path.exists(sample_corr_dir):
                os.mkdir(sample_corr_dir)
            sample_corr_df = perform_sample_corr(omics1, omics2, sample_corr_dir)
            # st.dataframe(sample_corr_df)

            ann_cls_cols = [i for i in data_type_dict.keys() if data_type_dict[i] in ['BIN','ORD']]
            ann_df2 = ann_df[ann_cls_cols]
            # st.dataframe(ann_df2)

            col_sample_corr_fig, col_sample_corr_settings = st.columns(2)
            with col_sample_corr_settings:
                st.write("Settings")
                meta_cols = st.multiselect("Select meta columns", ann_cls_cols)
            with col_sample_corr_fig:
                # st.dataframe(sample_corr_df)
                fig = plot_sample_corr(sample_corr_df, ann_df2, cols=meta_cols)
                st.pyplot(fig)

            gene_sample_corr_df = perform_calculate_gene_sample_corr(sample_corr_df, omics1, omics2, out_dir=sample_corr_dir)
            gene_sample_corr_df2 = perform_add_gene_name(gene_sample_corr_df, gene_map)
            
            gene_sample_corr1 = gsea_prerank(gene_sample_corr_df2, col_name="Corr_omics1", 
                out_dir=sample_corr_dir, processes=10,permutation_num=100)
            gene_sample_corr2 = gsea_prerank(gene_sample_corr_df2, col_name="Corr_omics2",
            out_dir=sample_corr_dir, processes=10,permutation_num=100)


            tab_table_sample_corr, tab_table_gene_sample_corr, tab_tabe_pathway_omics1, tab_tabe_pathway_omics2 = st.tabs(["Sample Correlation", "Gene-Sample Correlation", "Pathway Omics1", "Pathway Omics2"])
            with tab_table_sample_corr:
                st.dataframe(sample_corr_df)
            with tab_table_gene_sample_corr:
                st.dataframe(gene_sample_corr_df)
            with tab_tabe_pathway_omics1:
                st.dataframe(gene_sample_corr1.res2d)
            with tab_tabe_pathway_omics2:
                st.dataframe(gene_sample_corr2.res2d)
            

            

        with tab3:
            sample_cluster_dir = os.path.join(out_dir,"omicsx_sample_clustering")
            if not os.path.exists(sample_cluster_dir):
                os.mkdir(sample_cluster_dir)

            from utils.omics_sample_cluster import calculate_gene_high_and_low
            gene_high, gene_low = calculate_gene_high_and_low(omics1, omics2, out_dir=sample_cluster_dir)

            from utils.omics_sample_cluster import newclustermap2
            col_sample_cluster_fig, col_sample_cluster_settings = st.columns(2)
            with col_sample_cluster_settings:
                st.write("Settings")

            with col_sample_cluster_fig:
                tab_gene_high_omics1,tab_gene_high_omics2, tab_gene_low_omics1, \
                    tab_gene_low_omics2, tab_cluster_summary, tab_cluster_comparison = st.tabs(
                        ["Gene High Omics1", "Gene High Omics2", 
                         "Gene Low Omics1", "Gene Low Omics2",
                         "Cluster Summary", "Cluster Comparison"])
                with tab_gene_high_omics1:
                    mm1=omics1.loc[gene_high,:]
                    my_col_colors1, fig1 = newclustermap2(mm1, 
                        out_dir=sample_cluster_dir, name="gene_high_omics1")
                    st.pyplot(fig1)
                with tab_gene_high_omics2:
                    mm2=omics2.loc[gene_high,:]
                    my_col_colors2, fig2 = newclustermap2(mm2, 
                        out_dir=sample_cluster_dir, name="gene_high_omics2")
                    st.pyplot(fig2)
                with tab_gene_low_omics1:
                    mm3=omics1.loc[gene_low,:]
                    my_col_colors3, fig3 = newclustermap2(mm3, 
                        out_dir=sample_cluster_dir, name="gene_low_omics1")
                    st.pyplot(fig3)
                with tab_gene_low_omics2:
                    mm4=omics2.loc[gene_low,:]
                    my_col_colors4, fig4 = newclustermap2(mm4, 
                        out_dir=sample_cluster_dir, name="gene_low_omics2")
                    st.pyplot(fig4)

                with tab_cluster_summary:
                    types = perform_merge_col_colors(my_col_colors1, my_col_colors2, my_col_colors3, my_col_colors4)
                    types.replace({ 'g' : 'green', 'r' : 'red', 'c' : 'cyan','m':'magenta', 'y':'yellow','b':'blue','k':'black'}, inplace=True)
                    types.to_csv(os.path.join(sample_cluster_dir,"clustering.txt"),sep="\t")
                    fig = plot_merge_col_colors(types, out_dir=sample_cluster_dir)
                    st.pyplot(fig)
                
                with tab_cluster_comparison:
                    fig_cluster_comparison  , df_cluster_comparison = perform_plot_cluster_comparison(types, out_dir=sample_cluster_dir)
                    st.pyplot(fig_cluster_comparison)

            tab_table1, tab_table2, tab_table3, tab_table4, tab_table_cluster_comparison \
                = st.tabs(["Gene High Omics1", "Gene High Omics2", "Gene Low Omics1", "Gene Low Omics2", "Cluster Comparison"])
            with tab_table1:
                st.dataframe(my_col_colors1)
            with tab_table2:
                st.dataframe(my_col_colors2)
            with tab_table3:
                st.dataframe(my_col_colors3)
            with tab_table4:
                st.dataframe(my_col_colors4)
            with tab_table_cluster_comparison:
                st.dataframe(df_cluster_comparison)

        with tab4:

            col_group, col_diff1, col_diff2 = st.columns(3)
            options  = ["Pathology"] +  list(ann_df.columns.values)
            with col_group:
                group_select = st.selectbox("Select a group",options)

            group_diff_options = []
            if group_select == "Pathology" and "Pathology" in set(readme.columns.values):
                temp = get_matching_rows(readme, path_a_select)
                group_diff_options = sorted(list(set(temp['Pathology'])))

                with col_diff1:
                    group_diff_select1 = st.selectbox('A',group_diff_options)
                with col_diff2:
                    group_diff_select2 = st.selectbox('B',[i for i in group_diff_options if i != group_diff_select1])
                

                gp_fn_a_a = get_matching_path_with_pathology(readme, path_a_select, group_diff_select1)
                gp_fn_a_b = get_matching_path_with_pathology(readme, path_a_select, group_diff_select2)
                gp_path_a_a = os.path.join(data_dir, gp_fn_a_a)
                gp_path_a_b = os.path.join(data_dir, gp_fn_a_b)

                gp_fn_b_a = get_matching_path_with_pathology(readme, path_b_select, group_diff_select1)
                gp_fn_b_b = get_matching_path_with_pathology(readme, path_b_select, group_diff_select2)
                gp_path_b_a = os.path.join(data_dir, gp_fn_b_a)
                gp_path_b_b = os.path.join(data_dir, gp_fn_b_b)


                omics_a_a = pd.read_csv(gp_path_a_a,sep="\t",header=0,index_col=0)
                omics_a_b = pd.read_csv(gp_path_a_b,sep="\t",header=0,index_col=0)
                omics_b_a = pd.read_csv(gp_path_b_a,sep="\t",header=0,index_col=0)
                omics_b_b = pd.read_csv(gp_path_b_b,sep="\t",header=0,index_col=0)

               # First rename the columns in the original dataframes
                omics_a_a.columns = [i + f'.{group_diff_select1}' for i in omics_a_a.columns]
                omics_a_b.columns = [i + f'.{group_diff_select2}' for i in omics_a_b.columns]
                omics_b_a.columns = [i + f'.{group_diff_select1}' for i in omics_b_a.columns]
                omics_b_b.columns = [i + f'.{group_diff_select2}' for i in omics_b_b.columns]

                # Store the renamed column names for the statistical test
                samples_a_a = omics_a_a.columns.tolist()
                samples_a_b = omics_a_b.columns.tolist()
                samples_b_a = omics_b_a.columns.tolist()
                samples_b_b = omics_b_b.columns.tolist()

                # Merge DataFrames by index (no need for suffix since columns are already renamed)
                omics1_merged = omics_a_a.join(omics_a_b, how='outer')
                omics2_merged = omics_b_a.join(omics_b_b, how="outer")


                # calculate statistic test
                omics1_diff, omics2_diff = perform_two_diff_analysis(omics1_merged, omics2_merged, 
                                                                     samples_a_a, samples_a_b, 
                                                                     samples_b_a, samples_b_b)


                # Merge the two diff results on index (features)
                merged_diff = pd.merge(omics1_diff, omics2_diff, 
                                    left_index=True, right_index=True,
                                    suffixes=('_omics1', '_omics2'))


                merged_diff['Log2FoldChange'] = merged_diff['Log2FC(median)_omics1'] - merged_diff['Log2FC(median)_omics2']

                perform_add_gene_name(merged_diff, gene_map)
                perform_add_gene_name(omics1_diff, gene_map)
                perform_add_gene_name(omics2_diff, gene_map)

                gd_job_dir = os.path.join(out_dir,"group_diff","gsea_prerank")
                sstemp2 = gsea_prerank(merged_diff,col_name="Log2FoldChange",out_dir=gd_job_dir,processes=10,permutation_num=100)


                # Adjust layout to prevent label cutoff
                plt.tight_layout()

                col_group_diff_fig, col_group_diff_settings = st.columns(2)
                
                with col_group_diff_settings:

                    tab_group_diff_settings, tab_group_diff_details = st.tabs(["Settings", "Details"])
                    with tab_group_diff_settings:
                        st.write("Settings")
                        # Optional: Display correlation coefficient
                    with tab_group_diff_details:
                        correlation = merged_diff['Log2FC(median)_omics1'].corr(merged_diff['Log2FC(median)_omics2'])
                        st.write(f"Correlation coefficient between log2FC values: {correlation:.3f}")


                        st.write(gp_path_a_a)
                        st.write(gp_path_a_b)
                        st.write(gp_path_b_a)
                        st.write(gp_path_b_b)

                with col_group_diff_fig:
                    fig = plot_group_diff(merged_diff, class_a_select, class_b_select, group_diff_select1, group_diff_select2)
                    st.pyplot(fig)


                tab_table_pathway, tab_table_merged_diff, tab_table_diff = st.tabs(["Pathway", "Merged Diff", "Omics Diff"])
                with tab_table_pathway:
                    st.dataframe(sstemp2.res2d)
                with tab_table_merged_diff:
                    st.dataframe(merged_diff)
                with tab_table_diff:
                    col_diff1, col_diff2 = st.columns(2)
                    with col_diff1:
                        st.dataframe(omics1_diff)
                    with col_diff2:
                        st.dataframe(omics2_diff)

                            