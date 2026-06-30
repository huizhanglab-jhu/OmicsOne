import streamlit as st
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TABLEAU_COLORS
from utils.fasta import get_gene_map
from utils.data import detect_variable_type
from utils.params import read_readme
import numpy as np

import matplotlib.pyplot as plt

from matplotlib.patches import Patch

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.colors import ListedColormap
import gc

import io

@st.cache_data
def perform_get_gene_map(path):
    gene_map = get_gene_map(path)
    return gene_map

@st.cache_data
def perform_corr(method="spearman", cnv_data=None, rna_data=None, protein_data=None):
    from utils.data import fast_rowwise_correlation, fast_rowwise_spearman
    if method == "spearman":
        corr_cnv_rna = fast_rowwise_spearman(cnv_data, rna_data)
        corr_cnv_protein = fast_rowwise_spearman(cnv_data, protein_data)
        return corr_cnv_rna, corr_cnv_protein
    elif method == "pearson":
        corr_cnv_rna = fast_rowwise_correlation(cnv_data, rna_data)
        corr_cnv_protein = fast_rowwise_correlation(cnv_data, protein_data)
        return corr_cnv_rna, corr_cnv_protein

def mutations_app():
    
    
    data_dir = st.session_state.data_dir
    fasta_path = st.session_state.fasta_path
    
    gene_map = perform_get_gene_map(fasta_path)
    
    # st.write(f'gene_map: {gene_map}')
    # st.write(fasta_path)
    out_dir = st.session_state.out_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    job_dir = os.path.join(out_dir, 'genomics_mutations')
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
    else:
        # st.title("Quaity Control")
        st.title("OmicsOne Mutations Analysis")
        
        files = readme[(readme["Class"]=="Genomics") & (readme["Experiment.Method"].str.contains("Somatic_mutation")) & (
            readme["Quant.Method"] == "binary"
        )]["Path"].unique()
        
        meta_files = readme[(readme["Class"]=="Other_Meta")]["Path"].unique()
        
        maf_files = readme[(readme["Class"]=="Genomics") & (readme["Quant.Method"]=="maf") & (
            readme["Experiment.Method"]=="Somatic_mutation")]["Path"].unique()
        
        
        if len(files) == 0:
            st.write("No mutation files found in the data directory.")
            st.stop()
            
        if len(meta_files) == 0:
            st.write("No meta files found in the data directory.")
            st.stop()
            
        if len(maf_files) == 0:
            st.write("No maf files found in the data directory.")
            
            
        select_cols = st.columns([1,1])
        with select_cols[0]:
            file_select = st.selectbox('file', files, key='mutation_file')
            
            if len(maf_files) > 0:
                sel_maf = st.selectbox('maf_file', maf_files, key='mutation_maf_file')
            else:
                sel_maf = None
        
        data = pd.read_csv(os.path.join(data_dir, file_select), sep="\t", index_col=0)
        
        sample_cols = data.columns.tolist()
        
        data2 = data.copy(deep=True)
        data2['NEW.NUM_MUT']  = (data2 != 0).mean(axis=1)

        data2['NEW.Gene'] = data2.index.map(lambda x: gene_map.get(x.split(".")[0], dict(gene=x)).get('gene', x))

        
        if sel_maf is not None:
            maf = pd.read_csv(os.path.join(data_dir, sel_maf), sep="\t", index_col=0)
        else:
            maf = None
        
        
        
        with select_cols[1]:
            meta_select = st.selectbox('meta_file', meta_files, key='mutation_meta_file')
        
        meta = pd.read_csv(os.path.join(data_dir, meta_select), sep="\t", index_col=0,header=[0,1])
        
        
        
        with st.expander("Data"):
            
            if maf is not None:
                tabs = st.tabs(["Data","Meta","MAF"])
            else:
                tabs = st.tabs(["Data","Meta"])
            
            with tabs[0]:
                st.dataframe(data2)
            
            with tabs[1]:
                st.dataframe(meta)
            
            if maf is not None:
                with tabs[2]:
                    st.dataframe(maf)
                
        mutations_results_tabs = st.tabs(["Heatmap"])
        
        with mutations_results_tabs[0]:
            mutations_heatmap_cols = st.columns([1,1])
            
            with mutations_heatmap_cols[0]:
                st.subheader("Settings")
                items = [i[0] for i in meta.columns.tolist() if i[1] in ['ORD','BIN']]
                
                items_map = dict([(i,(i,j)) for i,j in meta.columns.tolist()])
                
                sel_items = st.multiselect("Select items to top annotations", items, default=items[:2])
                
                mut_ratio_threshold = st.number_input("Mutations Ratio Threshold", min_value=0.0, max_value=1.0, value=0.15, step=0.01)
                
                radio_maf = st.radio("Select MAF", ["No MAF", "MAF"], index=0, horizontal=True)
                
            data3 = data2.sort_values(by='NEW.NUM_MUT', ascending=False)
            data3 = data3[data3['NEW.NUM_MUT'] >= mut_ratio_threshold]
            data3 = data3.set_index('NEW.Gene')
            data3 = data3.loc[:, sample_cols]
            
            pairs = []
            for i in data3.columns:
                name = ''.join([f'{k}' for k in list(data3[i])])
                pairs.append((i,name))
                pairs = sorted(list(pairs),key=lambda x:x[1], reverse=True)
                
            df_sorted = data3[[i[0] for i in pairs]]
            df_sorted = df_sorted.map(float)
            
            df_sorted.index = [i.split(',')[0] for i in df_sorted.index]
            
            df_encoded = df_sorted.copy(deep=True)
            
            if radio_maf == "MAF":
                maf2 = maf[maf['Hugo_Symbol'].isin(df_sorted.index)]
                mut_types = {}
                for index,row in maf2.iterrows():
                    gene = row["Hugo_Symbol"]
                    var_cls  = row["Variant_Classification"]
                    sample = row["Tumor_Sample_Barcode"]
                    key = (sample, gene)
                    if key not in mut_types:
                        mut_types[key] = set()
                    mut_types[key].add(var_cls)
                
                df_var = df_sorted.copy(deep=True)
                for index,row in df_var.iterrows():
                    for sample in df_var.columns:
                        if row[sample] == 1:
                            key = (sample,index)
                            if key in mut_types:
                                var_cls = mut_types[key]
                                if len(var_cls) == 1:
                                    df_var.at[index,sample] = list(var_cls)[0]
                                else:
                                    df_var.at[index,sample] = "Multiple_Mutations"
                            else:
                                df_var.at[index,sample] = "Unknown"
                df_var.replace(0, "No_Mutation", inplace=True)
                
                unique_mutations = pd.unique(df_var.values.ravel())
                unique_mutations = sorted([x for x in unique_mutations if pd.notnull(x)])  # 去除 NaN
                
                # 为每种 mutation 分配颜色
                palette = sns.color_palette("Set2", len(unique_mutations))
                color_map = dict(zip(unique_mutations, palette))
                color_map["No_Mutation"] = "white"  # Gray for No_Mutation
                color_map["Frame_Shift_Del"] = "red"  # Red for Frame_Shift_Del
                default_color_map = {
                    "Frame_Shift_Del": "#d62728",           # red
                    "Frame_Shift_Ins": "#e377c2",           # pink
                    "In_Frame_Del": "#bcbd22",              # yellow-green
                    "In_Frame_Ins": "#7f7f7f",              # gray
                    "Missense_Mutation": "#17becf",         # cyan
                    "Nonsense_Mutation": "#ff9896",         # light red
                    "Splice_Site": "#9edae5",               # light cyan
                    "Translation_Start_Site": "#1f77b4",    # blue
                    "Multiple_Mutations": "#ff7f0e"         # orange
                }
                for i in default_color_map:
                    color_map[i] = default_color_map[i]

                
                number_map = {}
                for i, j in enumerate(unique_mutations):
                    number_map[j] = i
                
                df_encoded = df_var.map(lambda x: number_map.get(x, x))
                
    
                
                
                # st.dataframe(df_encoded)
                
            with mutations_heatmap_cols[1]:
                
                mut_heatmap_tabs = st.tabs(["Heatmap","Barplot"])
                # st.subheader("figure")
                with mut_heatmap_tabs[0]:
                # build col_colors
                
                    rows = []
                    
                    for sample in df_sorted.columns:
                        values = []
                        for i in sel_items:
                            col = items_map.get(i)
                            value = meta.loc[sample,col]
                            values.append(value)
                        rows.append(values)
                    col_df = pd.DataFrame(rows, columns=sel_items, index=df_sorted.columns)
                    col_df = col_df.replace(np.nan, "NA")
                    # st.dataframe(col_df)
                    
                    from utils.display import categorical_colormap
                    col_color_df = categorical_colormap(col_df)
                    
                    # st.dataframe(col_color_df)
                                    


                    # Define custom colors for 0 and 1
                    custom_colors = ["#d3d3d3", "#1f77b4"]  # Gray for 0, Blue for 1

                    # Create a custom colormap
                    
                    cmap = ListedColormap(custom_colors)
                    if radio_maf == "MAF":
                        ordered_colors = [color_map.get(i, "#d3d3d3") for i in unique_mutations]
                        cmap = ListedColormap(ordered_colors)
                    # Plot heatmap
                    # plt.figure(figsize=(20, 4))  # Set the figure size
                    
                    g = sns.clustermap(df_encoded, cmap=cmap, annot=False,
                                    col_cluster=False, row_cluster=False,
                                    col_colors=col_color_df, 
                                    figsize=(10,12),
                                        linewidths=0.5,       # Set grid line width
                        linecolor='grey',    # Set grid line color
                                    cbar_pos=None)  # Heatmap with annotations

                    # Access the `col_colors` and heatmap axes
                    col_colors_ax = g.ax_col_colors
                    heatmap_ax = g.ax_heatmap

                    # Adjust the `col_colors` axis height
                    col_colors_ax.set_position([col_colors_ax.get_position().x0,
                                                col_colors_ax.get_position().y0,
                                                col_colors_ax.get_position().width,
                                                col_colors_ax.get_position().height * 2])  # Double the height

                    # Adjust the heatmap size if needed
                    heatmap_ax.set_position([heatmap_ax.get_position().x0,
                                            heatmap_ax.get_position().y0 ,  # Shift down slightly
                                            heatmap_ax.get_position().width,
                                            heatmap_ax.get_position().height * 1])  # Slightly reduce height


                    # Add legend for column annotations
                    if radio_maf != "MAF":
                        legend_elements = [
                            Patch(facecolor="#1f77b4", label="True"),
                            Patch(facecolor="#d3d3d3", label="Fase")
                        ]
                        legend1 = plt.legend(
                            handles=legend_elements,
                            title="Mutation",
                            loc="upper left",
                            bbox_to_anchor=(1.15, 0.3),
                            bbox_transform=plt.gcf().transFigure
                        )
                        
                        plt.gca().add_artist(legend1)
                    else:
                        legend_elements = [
                            Patch(facecolor=color, label=label)
                            for label, color in color_map.items()
                        ]

                        # 添加 legend
                        legend =  plt.legend(
                            handles=legend_elements,
                            title="Mutation",
                            loc="upper left",
                            bbox_to_anchor=(1.15, 0.3),
                            bbox_transform=plt.gcf().transFigure
                        )
                        plt.gca().add_artist(legend)
        
                    
                    from utils.display import create_categorical_legend
                    for i, item in enumerate(sel_items):
                        legend = create_categorical_legend(col_df, 
                                                        col_color_df, 
                                                        item,
                                                        legend_title=item,
                                                        bbox_to_anchor=(1.15, 0.5 + i * 0.2))
                        plt.gca().add_artist(legend)


                    # Add both legends to the figure
                    # plt.gca().add_artist(legend3)  # Ensure the first legend is retained
                    # plt.gca().add_artist(legend2)
                    # Remove xticks labels
                    temp = g.ax_heatmap.set_xticks([])

                    # Move yticks to the left
                    g.ax_heatmap.yaxis.tick_left()
                    g.ax_heatmap.yaxis.set_label_position("left")
                    plt.subplots_adjust(right=1.0)
                    
                    # 3. 获取 heatmap 的位置，用于对齐右边 bar chart
                    heatmap_pos = heatmap_ax.get_position()
                    
                    # 4. 创建新的 axes 在 heatmap 右边
                    bar_width = 0.1  # 控制 bar chart 宽度
                    bar_ax = g.fig.add_axes([
                        heatmap_pos.x1 + 0.01,         # x 起点（紧贴 heatmap 右侧）
                        heatmap_pos.y0 ,                # y 起点
                        bar_width,                     # 宽度
                        heatmap_pos.height             # 高度匹配 heatmap
                    ])
                    
                    mutation_pct = (df_sorted.sum(axis=1) / df_sorted.shape[1]) * 100
                    

                    
                    # 5. 画横向 bar chart（基因在 y 轴，数值在 x 轴）
                    bar_ax.barh(
                        y=np.arange(len(df_sorted.index)),
                        width=mutation_pct.values,
                        color="skyblue",
                        edgecolor="black"
                    )
                    # Add labels manually
                    for i, v in enumerate(mutation_pct.values):
                        bar_ax.text(
                            v + 1,       # x-position: just past the end of the bar
                            i,           # y-position: bar’s center
                            f"{v:.1f}%", # label text
                            va="center"  # vertical alignment
                        )
                    
                    # 6. 对齐方向和格式
                    
                    bar_ax.set_ylim(-0.5,len(df_sorted.index)-0.5)  # y轴范围与 heatmap 匹配
                    bar_ax.invert_yaxis()                      # y轴方向匹配 heatmap
                    bar_ax.set_xlim(0, 100)
                    # 
                    bar_ax.set_xlabel("% samples\nwith mutations")
                    bar_ax.yaxis.set_visible(False)            # 不显示 y 轴标签
                    sns.despine(ax=bar_ax, left=True, right=True, top=True)


                    # out_path = os.path.join(job_dir, "mutation_hetamp.png")
                    # plt.savefig(out_path, dpi=300, bbox_inches='tight')
                    st.pyplot(g.figure,bbox_inches='tight')
                    
                    fig = g.figure
                    fig.set_size_inches(10, 8) 
                    
                    # fig.tight_layout()      # or g.figure.tight_layout()

                    # 2) write to a PDF buffer with a tight bbox
                    pdf_buf = io.BytesIO()
                    fig.savefig(pdf_buf,
                                format="pdf",
                                bbox_inches="tight",
                                pad_inches=0.1)   # optional padding around the edges
                    pdf_buf.seek(0)

                    # 两列布局，比例 9 : 1
                    col1, col2 = st.columns([3, 1])
                    with col2:
                        # 4. 添加下载按钮
                        st.download_button(
                            label="Download PDF",
                            data=pdf_buf,
                            file_name="plot.pdf",
                            mime="application/pdf"
                        )

                with mut_heatmap_tabs[1]:
                    if radio_maf != "MAF":
                        st.write("Barplot is not available when MAF is not selected.")
                        st.stop()
                    rows = []
                    from collections import Counter
                 
                    target_mutations = [i for i in unique_mutations if i != "No_Mutation"]
                    for index,row in df_var.iterrows():
                        counts = Counter(list(row))
                        new_row = [counts.get(i,0) for  i in target_mutations]
                        rows.append(new_row)
                    df_counts = pd.DataFrame(rows, columns=target_mutations, index=df_var.index)
                    
                    # Convert counts to percentages row-wise
                    df_percent = df_counts.div(df_counts.sum(axis=1), axis=0) * 100
                    
                    # Plot stacked horizontal bar
                    fig, ax = plt.subplots(figsize=(6, 8),dpi=300)

                    left = pd.Series([0] * df_percent.shape[0], index=df_percent.index)
                    for mut in df_percent.columns:
                        ax.barh(df_percent.index, df_percent[mut], left=left, color=color_map.get(mut, "#ccc"), label=mut)
                        left += df_percent[mut]

                    # Style
                    ax.set_xlabel("Percentage (%)")
                    ax.set_xlim(0, 100)
                    ax.set_ylim(-0.5, len(df_percent.index) - 0.5)
                    ax.invert_yaxis()  # optional, to match heatmap ordering
                    ax.legend(title="Mutation Type", bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.tight_layout()
                    
                    st.pyplot(fig)
                    
                    from utils.display import create_figure_pdf_buf
                    # 两列布局，比例 9 : 1
                    col1, col2 = st.columns([3, 1])
                    with col2:
                        # 4. 添加下载按钮
                        st.download_button(
                            label="Download PDF",
                            data=create_figure_pdf_buf(fig),
                            file_name="plot.pdf",
                            mime="application/pdf",
                            key="download_mut_percent_barplot_pdf"
                        )

                
                
    
    
def cnv_app():
    st.title("OmicsOne CNV Analysis")
    data_dir = st.session_state.data_dir
    fasta_path = st.session_state.fasta_path
    
    gene_map = perform_get_gene_map(fasta_path)
    st.write(f"Using fasta file: {fasta_path}")
    
    chrom_path = st.session_state.chrom_path
    st.write(f"Using chromosome file: {chrom_path}")
    
    chrom = pd.read_csv(chrom_path, sep="\t", header=0)
    
    
    n = 0

    rows = []
    for index,row in chrom.iterrows():
        row['Start'] = n
        n += int(str(row["Total length (bp)"]).replace(",",""))
        rows.append(row)
        
    chrom = pd.DataFrame(rows)
    chrom['Chromosome'] = chrom['Chromosome'].astype(str)
    
    total_bp =  np.sum([int(str(i).replace(",","")) for i in chrom['Total length (bp)']])
    # st.write(total_bp)
    starts = chrom['Start'].tolist()
    
    chr_bounds = starts + [int(total_bp)]
    
    chr_names = ["\n"+str(i) if (i % 2== 0) and (i> 10) else str(i) for i in range(1,23)] + ["X","Y"]
    chr_names2 = [str(i) + "    " if (i % 2== 0) and (i> 10) else str(i) for i in range(1,23)] + ["X","Y"]
    
    chr_ticks = [int((chr_bounds[i] + chr_bounds[i+1])/2) for i in range(len(chr_bounds)-1)]
    
    st.write(f"Total length of all chromosomes: {total_bp}")
    
    start_map = dict(zip(chrom["Chromosome"],chrom["Start"]))
    
    cytoband_path = st.session_state.cytoband_path
    
    from utils.genome import get_cytoband_map
    cytoband_d, cytobands = get_cytoband_map(cytoband_path)
    
    from utils.genome import add_band_info
    gene_map = add_band_info(gene_map, cytobands)
    
    # for i in start_map:
    #     st.write(i)
    #     st.write(start_map[i])
    #     break
    
    # for i in gene_map:
    #     st.write(i)
    #     st.write(gene_map[i])
    #     break
    
    def get_location(gene_id):
        try:
            d = gene_map[str(gene_id).split(".")[0]]
            chr = d['chr']
            pos = d['offset']
            return start_map[chr] + int(pos)
        except:
            return np.nan
        
    def get_gene(gene_id):
        try:
            d = gene_map[str(gene_id).split(".")[0]]
            return d['gene']
        except:
            return np.nan
        
    def get_band(gene_id):
        try:
            d = gene_map[str(gene_id).split(".")[0]]
            return d['band']
        except:
            return np.nan
        
    def get_chr(gene_id):
        try:
            d = gene_map[str(gene_id).split(".")[0]]
            return d['chr']
        except:
            return np.nan
    
    # select protein file
    readme = read_readme(data_dir)
    if readme is None:
        st.write(f'Fail to read in readme.xlsx from {data_dir}')
        st.stop()
    
    cnv_file_cols = st.columns([1,1])    
    with cnv_file_cols[0]:
        cnv_files = readme[(readme["Class"]=="Genomics") & 
                           (readme["Experiment.Method"]=="WES_CNV") & 
                           (readme["Quant.Method"]=="ratio")]["Path"].unique()
        sel_cnv_file = st.selectbox("Select CNV file (log2ratio)", cnv_files, key="cnv_file")
        
    with cnv_file_cols[1]:
        gistic_files = readme[(readme["Class"]=="Genomics") & 
                           (readme["Experiment.Method"]=="WES_CNV") & 
                           (readme["Quant.Method"]=="gistic")]["Path"].unique()
        sel_gistic_file = st.selectbox("Select CNV file (gistic)", gistic_files, key="gistic_file")
    
    other_file_cols = st.columns([1,1])
    with other_file_cols[0]:
        condition = (
        ((readme["Class"] == "Transcriptomics") & (readme["Experiment.Method"] == "RNASeq") & (readme["Group"] == "gene") ) 
        )

        rna_files = readme[condition]["Path"].unique()
        sel_rna_file = st.selectbox("Select RNA file", rna_files, key="rna_file")
    
    with other_file_cols[1]:
        condition = (
        ((readme["Class"] == "Proteomics") & (readme["Experiment.Method"] == "MS") & (readme["Group"] == "gene") ) 
        )
        
        protein_files = readme[condition]["Path"].unique()
        sel_protein_file = st.selectbox("Select Protein file", protein_files, key="protein_file")
    
    cnv_data = pd.read_csv(os.path.join(data_dir, sel_cnv_file), sep="\t", index_col=0)
    gistic_data = pd.read_csv(os.path.join(data_dir, sel_gistic_file), sep="\t", index_col=0)
    rna_data = pd.read_csv(os.path.join(data_dir, sel_rna_file), sep="\t", index_col=0)
    protein_data = pd.read_csv(os.path.join(data_dir, sel_protein_file), sep="\t", index_col=0)
    
    with st.expander("Data"):
        tabs = st.tabs(["CNV(log2ratio)","CNV(gistic)","RNA","Protein","Chromosome"])
        with tabs[0]:
            st.dataframe(cnv_data)
        with tabs[1]:
            st.dataframe(gistic_data)
        with tabs[2]:
            st.dataframe(rna_data)
        with tabs[3]:
            st.dataframe(protein_data)
        with tabs[4]:
            st.dataframe(chrom)
            
    with st.expander("Settings"):
        settings_tabs = st.tabs(["Preprocessing","Figure"])
        with settings_tabs[0]:
            st.subheader("Preprocessing")
            threshold_cols = st.columns([1,1,1])
            cnv_threshold = 0.0
            rna_threshold = 0.0
            protein_threshold = 0.0
            # with threshold_cols[0]:
            #     cnv_threshold = st.number_input(
            #         "CNV Threshold (0-1), 0: no missing values, 1: no requiremet", min_value=0.0, max_value=1.0, value=float(0), step=0.01)
            # with threshold_cols[1]:
            #     rna_threshold = st.number_input(
            #         "RNA Threshold (0-1), 0: no missing values, 1: no requiremet", min_value=0.0, max_value=1.0, value=float(0), step=0.01)
            # with threshold_cols[2]:
            #     protein_threshold = st.number_input(
            #         "Protein Threshold (0-1), 0: no missing values, 1: no requiremet", min_value=0.0, max_value=1.0, value=float(0), step=0.01)
            with threshold_cols[0]:
                corr_value_threshold = st.number_input(
                "Correlation Threshold", min_value=0.0, max_value=1.0, value=float(0.5), step=0.01)
    from utils.data import replace_inf_to_nan, remove_rows_by_missing_values, convert_to_float
    
    cnv_data = cnv_data.map(float)
    rna_data = rna_data.map(float)
    protein_data = protein_data.map(float)
    
    cnv_data_cleaned = (
        cnv_data
        .pipe(convert_to_float)
        .pipe(replace_inf_to_nan)
        .pipe(remove_rows_by_missing_values, threshold=cnv_threshold) 
    )
    
    rna_data_cleaned = (
        rna_data
        .pipe(convert_to_float)
        .pipe(replace_inf_to_nan)
        .pipe(remove_rows_by_missing_values, threshold=rna_threshold) 
    )
    
    protein_data_cleaned = (
        protein_data
        .pipe(convert_to_float)
        .pipe(replace_inf_to_nan)
        .pipe(remove_rows_by_missing_values, threshold=protein_threshold) 
    )
    
    common_genes = set(cnv_data_cleaned.index) & set(rna_data_cleaned.index) & set(protein_data_cleaned.index)
    common_samples = set(cnv_data_cleaned.columns) & set(rna_data_cleaned.columns) & set(protein_data_cleaned.columns)
    
    common_genes = sorted(list(common_genes))
    common_samples = sorted(list(common_samples))
    
    cnv_for_corr = cnv_data_cleaned.loc[common_genes, common_samples]
    rna_for_corr = rna_data_cleaned.loc[common_genes, common_samples]
    protein_for_corr = protein_data_cleaned.loc[common_genes, common_samples]
    
    # from utils.data import fast_rowwise_correlation, fast_rowwise_spearman
    # corr_cnv_rna = fast_rowwise_spearman(cnv_for_corr, rna_for_corr)
    # corr_cnv_protein = fast_rowwise_spearman(cnv_for_corr, protein_for_corr)
    
    corr_cnv_rna, corr_cnv_protein = \
        perform_corr(method="spearman",
                     cnv_data=cnv_for_corr, 
                     rna_data=rna_for_corr, 
                     protein_data=protein_for_corr)
    
    corr_cnv_rna["CNV"] = corr_cnv_rna.index
    corr_cnv_rna_long = pd.melt(corr_cnv_rna, id_vars=["CNV"], 
                                var_name="RNA", value_name="Correlation")

    corr_cnv_rna_long_filtered = \
        corr_cnv_rna_long[abs(corr_cnv_rna_long["Correlation"]) > corr_value_threshold]
        
    # del corr_cnv_rna_long
    # gc.collect()
    
    def build_band(row):
        chrom = row['cnv.chr']
        band = row['cnv.cytoband']
        key = f'{chrom}{band[0]}'
        return key
    
    corr_cnv_rna_long_filtered['cnv.location'] = corr_cnv_rna_long_filtered['CNV'].map(get_location)
    corr_cnv_rna_long_filtered['rna.location'] = corr_cnv_rna_long_filtered['RNA'].map(get_location)
    corr_cnv_rna_long_filtered['cnv.gene'] = corr_cnv_rna_long_filtered['CNV'].map(get_gene)
    corr_cnv_rna_long_filtered['rna.gene'] = corr_cnv_rna_long_filtered['RNA'].map(get_gene)
    
    corr_cnv_rna_long_filtered["cnv.cytoband"] = corr_cnv_rna_long_filtered['CNV'].map(get_band)
    corr_cnv_rna_long_filtered["cnv.chr"] = corr_cnv_rna_long_filtered['CNV'].map(get_chr)

    corr_cnv_rna_long_filtered["rna.cytoband"] = corr_cnv_rna_long_filtered['RNA'].map(get_band)
    corr_cnv_rna_long_filtered["rna.chr"] = corr_cnv_rna_long_filtered['RNA'].map(get_chr)
    
    # st.write(corr_cnv_rna_long_filtered.shape)
    corr_cnv_rna_long_filtered = corr_cnv_rna_long_filtered.dropna()
    # st.write(corr_cnv_rna_long_filtered.shape)
    
    corr_cnv_rna_long_filtered['cnv.band'] = corr_cnv_rna_long_filtered.apply(build_band,axis=1)
    
    corr_cnv_protein["CNV"] = corr_cnv_protein.index
    corr_cnv_protein_long = pd.melt(corr_cnv_protein, id_vars=["CNV"], 
                                var_name="Protein", value_name="Correlation")
    
    corr_cnv_protein_long_filtered = \
        corr_cnv_protein_long[abs(corr_cnv_protein_long["Correlation"]) > corr_value_threshold]
        
    # del corr_cnv_protein_long
    # gc.collect()
    
    corr_cnv_protein_long_filtered['cnv.location'] = corr_cnv_protein_long_filtered['CNV'].map(get_location)
    corr_cnv_protein_long_filtered['protein.location'] = corr_cnv_protein_long_filtered['Protein'].map(get_location)
    corr_cnv_protein_long_filtered['cnv.gene'] = corr_cnv_protein_long_filtered['CNV'].map(get_gene)
    corr_cnv_protein_long_filtered['protein.gene'] = corr_cnv_protein_long_filtered['Protein'].map(get_gene)
    
    corr_cnv_protein_long_filtered["cnv.cytoband"] = corr_cnv_protein_long_filtered['CNV'].map(get_band)
    corr_cnv_protein_long_filtered["cnv.chr"] = corr_cnv_protein_long_filtered['CNV'].map(get_chr)

    corr_cnv_protein_long_filtered["protein.cytoband"] = corr_cnv_protein_long_filtered['Protein'].map(get_band)
    corr_cnv_protein_long_filtered["protein.chr"] = corr_cnv_protein_long_filtered['Protein'].map(get_chr)
    
    corr_cnv_protein_long_filtered = corr_cnv_protein_long_filtered.dropna()
    
    corr_cnv_protein_long_filtered['cnv.band'] = corr_cnv_protein_long_filtered.apply(build_band,axis=1)
        
    
    # for i in gene_map:
    #     st.write(i)
    #     st.write(gene_map[i])
    #     break
    
    with st.expander("Results"):
        result_tabs = st.tabs(
            ["Common Genes","Common Samples","corr_cnv_rna","corr_cnv_protein"])
        with result_tabs[0]:
            st.subheader(f"Common Genes: {len(common_genes)}")
            st.dataframe(pd.DataFrame(sorted(list(common_genes)), columns=["Gene"]))
        with result_tabs[1]:
            st.subheader(f"Common Samples: {len(common_samples)}")
            st.dataframe(pd.DataFrame(sorted(list(common_samples)), columns=["Sample"]))
        with result_tabs[2]:
            st.subheader(f"Correlation between CNV and RNA")
            st.dataframe(corr_cnv_rna)
            st.write(corr_cnv_rna_long_filtered.shape)
            st.dataframe(corr_cnv_rna_long_filtered.head(10))
            temp = corr_cnv_rna_long
            st.write(temp[(temp["CNV"]=="ENSG00000000003.15")&(temp["RNA"].str.contains("00000126953"))])
        with result_tabs[3]:
            st.subheader(f"Correlation between CNV and Protein")
            st.dataframe(corr_cnv_protein)
            st.write(corr_cnv_protein_long_filtered.shape)
            st.dataframe(corr_cnv_protein_long_filtered.head(10))

        
    with st.expander("Figure"):
        figure_cols = st.columns([1,1])
        from plots.cnv import draw_corr
        with figure_cols[0]:
            # st.dataframe(corr_cnv_rna_long_filtered.head(10))
            # st.write(chr_bounds)
            # st.write(chr_names)
            # st.write(chr_ticks)
            # st.write(chr_names2)
            # st.write(start_map)
            fig = draw_corr(corr_cnv_rna_long_filtered, 
                            genes = common_genes, 
                            samples = common_samples,
                            chr_bounds = chr_bounds,
                            chr_names = chr_names,
                            chr_ticks = chr_ticks,
                            chr_names2 = chr_names2,
                            total_bp = total_bp,
                            gene_map = gene_map,
                            start_map = start_map,
                            cytoband_d = cytoband_d,
                            gistic_data = gistic_data,
                            corr_with = "RNA", 
                      )
            st.pyplot(fig)
            
        with figure_cols[1]:
            fig = draw_corr(corr_cnv_protein_long_filtered, 
                            genes = common_genes, 
                            samples = common_samples,
                            chr_bounds = chr_bounds,
                            chr_names = chr_names,
                            chr_ticks = chr_ticks,
                            chr_names2 = chr_names2,
                            total_bp = total_bp,
                            gene_map = gene_map,
                            start_map = start_map,
                            cytoband_d = cytoband_d,
                            gistic_data = gistic_data,
                            corr_with = "Protein",
                      )
            st.pyplot(fig)
                
    

    
    