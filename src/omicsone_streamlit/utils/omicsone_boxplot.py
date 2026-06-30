import streamlit as st
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns   
from scipy import stats

def boxplot_page(readme: pd.DataFrame,gene_map: dict,data_dir: str):

    class_options = [i for i in set(readme['Class']) if not re.search('Other.',i)]
    class_select = st.selectbox('class',class_options)       

    col_a, col_b = st.columns(2)
    
    group_a_options = list(readme[readme['Class']==class_select]['Path'])

    with col_a:
        group_a_select = st.selectbox('A',group_a_options)         

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
        group_b_select = st.selectbox('B',group_b_options)
    row_b = readme[readme['Path']==group_b_select].iloc[0]
    b_Pathology = row_b["Pathology"]

    file_a  = os.path.join(data_dir, group_a_select)
    a_df = pd.read_csv(file_a,sep="\t",index_col=0)
    # st.dataframe(a_df)

    file_b = os.path.join(data_dir, group_b_select)
    b_df = pd.read_csv(file_b,sep="\t",index_col=0)
    # st.dataframe(b_df)

    entries = sorted(list(set(a_df.index) & set(b_df.index)))
    

    # according to the class_select, get the gene name
    if class_select in ['Phosphoproteomics']:
        entries = [i.split('|')[0] for i in entries]
        entries = [f"{i}, {gene_map.get(i.split('.')[0])['gene'] if i.split('.')[0] in gene_map else 'NA'}" for i in entries]
    elif class_select in ['Transcriptomics'] and entries[0].split('_')[0] == 'circ':
        entries = [f"{i},{gene_map.get(i.split('_')[-1].split('.')[0])['gene'] if i.split('_')[-1].split('.')[0] in gene_map else 'NA'}" for i in entries]
    else:
        entries = [f"{i},{gene_map.get(i.split('.')[0])['gene'] if i.split('.')[0] in gene_map else 'NA'}" for i in entries]

    gene_select = st.selectbox('Gene',entries)
    gene_select, gene_select_symbol = gene_select.split(',')

    rows = []
    samples_a = list(a_df.columns.values)

    row_a = a_df[a_df.index==gene_select].iloc[0]

    for i in samples_a:
        rows.append([i,row_a[i],a_Pathology])

    samples_b = list(b_df.columns.values)

    row_b = b_df[b_df.index==gene_select].iloc[0]

    for i in samples_b:
        rows.append([i,row_b[i], b_Pathology])

    data = pd.DataFrame(rows,columns=['Sample','Value','Group'])
    # st.dataframe(data)

    group_values = [data[data['Group'] == group]['Value'] for group in [a_Pathology, b_Pathology]]
    # st.write(a_Pathology, b_Pathology)
    # st.write(group_values)
    stat, pvalue = stats.mannwhitneyu([i for i in group_values[0] if not pd.isna(i)], 
                                        [ i for i in group_values[1] if not pd.isna(i)],
                                        alternative='two-sided')

    st.write(f"p-value: {pvalue}")

    fig,ax = plt.subplots(figsize=(3,4))
    
    sns.boxplot(x="Group", y="Value", data=data,hue='Group')
    # sns.swarmplot(x="Group", y="Value", data=data, color=".25", hu)

    # Add significance annotation
    y_max = data['Value'].max()
    y_min = data['Value'].min()
    y_range = y_max - y_min
    
    # Add significance bar and star
    if pvalue < 0.001:
        sig_symbol = '***'
    elif pvalue < 0.01:
        sig_symbol = '**'
    elif pvalue < 0.05:
        sig_symbol = '*'
    else:
        sig_symbol = 'ns'

    ax.plot([0, 0, 1, 1], [y_max + y_range*0.1, y_max + y_range*0.15, 
            y_max + y_range*0.15, y_max + y_range*0.1], 'k-', linewidth=1)
    ax.text(0.5, y_max + y_range*0.13, sig_symbol, ha='center', va='bottom')


    ax.set_title(f"{gene_select_symbol}")

    ax.set_ylabel('Log2 abundance')
    ax.set_xlabel('')

    col1, col2 = st.columns(2)
    with col1:
        st.pyplot(fig)

