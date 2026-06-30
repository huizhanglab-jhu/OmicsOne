import os,sys,re
import pandas as pd
import numpy as np
import streamlit as st
from tqdm import tqdm

import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


@st.cache_data
def calculate_missing_values(data_df):
    rows = []
    for index,row in tqdm(data_df.iterrows()):
        missing = row.isnull().sum()
        missing_ratio = missing / len(row)
        averge_intensity = np.mean([i for i in row if not np.isnan(i)])
        rows.append([averge_intensity, missing, missing_ratio])
    result = pd.DataFrame(rows, columns=['Average Intensity', 'Missing Values', 'Missing Rate'])
    return result
                

@st.cache_data
def plot_missing_values(df):
    cols = [i for i in df.columns.values if i != "Intensity.Reference"]
    data_df = df[cols]
    data_df = data_df.replace('None',np.nan) \
        .map(float) \
        
    result = calculate_missing_values(data_df)
    # Create the KDE plot using fig, ax
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.kdeplot(
        data=result, 
        x="Average Intensity", 
        y="Missing Rate", 
        fill=True, 
        cmap="Blues",
        thresh=0,
        levels=100,
        ax=ax
    )
    ax.set_title('Missing Mechanism')
    ax.set_xlabel('Mean Abundance')
    ax.set_ylabel('Missing Rate')
    
    st.pyplot(fig)
    return fig

@st.cache_data
def build_uniq(df, meta_df, ):
    # read sample file
    uniq_sets = meta_df['Set'].unique()
    
    sample_set_map = dict([(index,row['Set']) for index,row in meta_df.iterrows()])

    samples = [index for index,row in meta_df.iterrows()]

    df = df.replace('None',np.nan) \
            .map(float) \
            .dropna(how='all')
    # print(df.shape)
    # for i in df.columns.values:
    #     print("!!!" + i + "!!!")
    
    set_seq_map = dict([(i,set()) for i in uniq_sets])

    for sample in tqdm(sample_set_map):
        peptides = list(df[pd.notna(df[sample])].index)
        sequences = [seq for site,gene,seq in peptides]
        sample_set = sample_set_map[sample]
        set_seq_map[sample_set].update(set(sequences))
    
        
    rows = []
    pep50 = []
    set_size = len(set_seq_map)
    for sample_set in tqdm(sorted(set_seq_map)):
        pep = set_seq_map[sample_set]
        uniqs = 0
        
        common_50 = 0
        common_100 = 0
        not_uniq = 0
        
        for p in pep:
            shares = 0
            uniq = True
            for s in sorted(set_seq_map):
                if s != sample_set and p in set_seq_map[s]:
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
    # print(len(set(pep50)))
    # print(result)
    return result
    
