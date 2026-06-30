#!/tomcat/python3env/bin/python3
import os,sys,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
from matplotlib.pyplot import figure, show
from scipy.stats import spearmanr, variation
import io
# import cgi
import base64
from collections import OrderedDict
import gseapy
import statsmodels.stats.multitest as smt

from tqdm import tqdm
import streamlit as st

def orig_loop_corr(omics1a, omics2a, omics1aflat, omics2aflat):
    corr  = []
    for index, row in tqdm(omics1a.iterrows()):
        corr1,p=spearmanr(omics1a.loc[index,:].T,omics2a.loc[index,:].T)
        if omics1aflat.quantile(0.1)[0]<0 and omics1aflat.quantile(1)[0]<200:
            cv1 = variation(np.power(2,omics1a.loc[index,:]), axis=0)
        else:
            cv1 = variation(omics1a.loc[index,:], axis=0)
        if omics2aflat.quantile(0.1)[0]<0 and omics2aflat.quantile(1)[0]<200:
            cv2 = variation(np.power(2,omics2a.loc[index,:]), axis=0)
        else:
            cv2 = variation(omics2a.loc[index,:], axis=0)
        corr.append([index,corr1,p,cv1,cv2])
    return corr

def loop_corr(omics1aflat, omics2aflat, omics1a, omics2a):
    from scipy.stats import spearmanr, variation
    import numpy as np
    from tqdm import tqdm

    # Precompute quantiles
    omics1a_q10 = omics1aflat.quantile(0.1)[0]
    omics1a_q100 = omics1aflat.quantile(1)[0]
    omics2a_q10 = omics2aflat.quantile(0.1)[0]
    omics2a_q100 = omics2aflat.quantile(1)[0]

    # Convert DataFrames to dictionaries for faster row access
    omics1a_dict = omics1a.to_dict(orient='index')
    omics2a_dict = omics2a.to_dict(orient='index')

    # Initialize result list
    corr = []

    # Iterate over indices and compute results
    for index in tqdm(omics1a_dict.keys()):
        # Access rows directly from dictionaries
        omics1a_row = np.array(list(omics1a_dict[index].values()))
        omics2a_row = np.array(list(omics2a_dict[index].values()))

        # Compute Spearman correlation
        corr1, p = spearmanr(omics1a_row, omics2a_row)

        # Compute coefficient of variation for omics1a
        if omics1a_q10 < 0 and omics1a_q100 < 200:
            cv1 = variation(np.power(2, omics1a_row), axis=0)
        else:
            cv1 = variation(omics1a_row, axis=0)

        # Compute coefficient of variation for omics2a
        if omics2a_q10 < 0 and omics2a_q100 < 200:
            cv2 = variation(np.power(2, omics2a_row), axis=0)
        else:
            cv2 = variation(omics2a_row, axis=0)

        # Append results to the list
        corr.append([index, corr1, p, cv1, cv2])
    return corr



def calculate_corr(path1, path2, out_dir, gene_map=None):

    omics1=pd.read_csv(path1,sep="\t",header=0,index_col=0)
    omics2=pd.read_csv(path2,sep="\t",header=0,index_col=0)
    omics1=omics1.apply (pd.to_numeric, errors='coerce').dropna()
    omics2=omics2.apply (pd.to_numeric, errors='coerce').dropna()
    omics1 = omics1.loc[~omics1.index.duplicated(keep='first'),~omics1.columns.duplicated(keep='first')]
    omics2 = omics2.loc[~omics2.index.duplicated(keep='first'),~omics2.columns.duplicated(keep='first')]
    # annotation=pd.read_csv(wo+os.sep+"omics3.txt",sep="\t",header=0,index_col=0)
    geneinorder=set(omics1.index).intersection(set(omics2.index))
    sampleinorder=set(omics1.columns).intersection(set(omics2.columns))


    if len(geneinorder)>9 and len(sampleinorder)>9:
        omics2a=omics2.loc[list(geneinorder),list(sampleinorder)]
        omics1a=omics1.loc[list(geneinorder),list(sampleinorder)]
        omics1a.to_csv(os.path.join(out_dir,"omics1a.txt"),sep="\t")
        omics2a.to_csv(os.path.join(out_dir,"omics2a.txt"),sep="\t")    
        # omics2az=pd.DataFrame(np.power(2,omics2a),columns=omics2a.columns, index=omics2a.index)
        omics2az=pd.DataFrame(omics2a,columns=omics2a.columns, index=omics2a.index)
        omics1az=pd.DataFrame(omics1a,columns=omics1a.columns, index=omics1a.index)

        omics1aflat=pd.DataFrame(omics1a.values.flatten())
        omics2aflat=pd.DataFrame(omics2a.values.flatten())
        corr= loop_corr(omics1aflat, omics2aflat, omics1a, omics2a)

        #     print([index,corr1,p])
        df_corrgenewise=pd.DataFrame(np.asarray(corr))
        df_corrgenewise.apply(pd.to_numeric, errors='ignore', downcast='float')
        df_corrgenewise.columns=["Gene","Gene Correlation","P","CV1","CV2"]
        df_corrgenewise['P'] = pd.to_numeric(df_corrgenewise['P'], errors='coerce')
        df_corrgenewise = df_corrgenewise[pd.notna(df_corrgenewise['P'])]
        df_corrgenewise["Gene Correlation"]=df_corrgenewise["Gene Correlation"].astype(float)
        df_corrgenewise.set_index('Gene', inplace=True)
        df_corrgenewise['BH adjusted P'] =  smt.multipletests(list(df_corrgenewise['P'].map(float)), method='fdr_bh')[1]
        df_corrgenewise.sort_values(by="Gene Correlation",inplace=True,ascending=False)
        if gene_map is not None:
            # print(gene_map)
            # print(df_corrgenewise.head(5))
            df_corrgenewise['Gene'] = [gene_map.get(i.split(".")[0],
                                                    {'gene':'NA'})['gene'] for i in df_corrgenewise.index]
        df_corrgenewise.to_csv(os.path.join(out_dir,"gene_wise_corr.txt"),sep="\t")
    return df_corrgenewise                           

def plot_corr(df_corrgenewise,out_dir):
    fig = plt.figure(figsize=(4,4))
    # Create the figure and axis object
    fig, ax = plt.subplots()
    sns.set_style("white")
    # ax=sns.distplot(df_corrgenewise["Gene Correlation"],color="orangered")
    # sns.displot(df_corrgenewise, x="Gene Correlation", color="orangered", kde=True)
    sns.histplot(df_corrgenewise["Gene Correlation"], color="orangered", stat='density',kde=True, ax=ax)
    plt.title('Gene-wise correlation')
    ax.set_xlabel("Spearman's correlation",size=11,rotation=0)
    ax.set_ylabel("Probability density",size=11,rotation=90)
    mediancorr=df_corrgenewise["Gene Correlation"].median()
    plt.text(-0.45, 1.4, "Median={0:.2f}".format(mediancorr),fontsize=12)
    plt.xticks(size=11)
    plt.yticks(size=11)
    plt.xlim((-0.5,1))
    # plt.ylim((0,2.3))
    plt.savefig(os.path.join(out_dir,"genecorrelation.png"),dpi=100,bbox_inches = 'tight')
    return fig

def gsea_prerank(df,col_name, out_dir, gene_sets='MSigDB_Hallmark_2020',processes=4,permutation_num=100):
    if 'Gene' in set(df.columns.values):
        df = df.drop_duplicates(subset='Gene')
        df = df.set_index('Gene')
    # print(df)
    sstemp2 = gseapy.prerank(rnk=df.loc[:,col_name], 
                         gene_sets=gene_sets,
                         processes=processes,
                         permutation_num=permutation_num, # reduce number to speed up test
                         outdir=out_dir,format='png',
                             )
    return sstemp2



