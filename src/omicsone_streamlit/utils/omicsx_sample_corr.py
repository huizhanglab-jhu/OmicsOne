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



def calculate_sample_corr(omics1, omics2, out_dir):

    corr=[]  
    for each in omics1.columns:
        corr1,p=spearmanr(omics1[each],omics2[each])
        corr.append([each,corr1,p])
    df_samplecorr=pd.DataFrame(np.asarray(corr))  
    df_samplecorr.columns=["Sample","Corr","p-value"]
    df_samplecorr.set_index('Sample', inplace=True)
    df_samplecorr["Corr"]=df_samplecorr["Corr"].astype(float)
    df_samplecorr.sort_values(by="Corr",inplace=True,ascending=False)
    df_samplecorr.to_csv(os.path.join(out_dir,"sample_wise_correlation.txt"),sep="\t")
    return df_samplecorr

def calculate_gene_sample_corr(df_samplecorr, omics1a, omics2a, out_dir):
    # correlation between gene expression and sample correlation
    corr=[]
    for index, row in omics2a.iterrows():
        corr1,p1=spearmanr(df_samplecorr.loc[omics2a.columns,"Corr"],omics1a.loc[index,:].T)
        corr2,p2=spearmanr(df_samplecorr.loc[omics2a.columns,"Corr"],omics2a.loc[index,:].T)
        corr.append([index,corr1,p1,corr2,p2])
    df_corrRNAproteinGene=pd.DataFrame(np.asarray(corr)) 
    df_corrRNAproteinGene.columns=["Gene","Corr_omics1","P_omics1","Corr_omics2","P_omics2"]
    df_corrRNAproteinGene.set_index('Gene', inplace=True)
    df_corrRNAproteinGene["Corr_omics1"]=df_corrRNAproteinGene["Corr_omics1"].astype(float)
    df_corrRNAproteinGene["Corr_omics2"]=df_corrRNAproteinGene["Corr_omics2"].astype(float)
    df_corrRNAproteinGene.sort_values(by="Corr_omics1",inplace=True,ascending=False)
    if os.path.exists(out_dir):
        df_corrRNAproteinGene.to_csv(os.path.join(out_dir,"sample_wise_corr.txt"),sep="\t")
    return df_corrRNAproteinGene
