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
from .diff import compare_two_groups

def get_matching_rows(readme, path_select):
    # Get parameters from selected path
    params = readme[readme['Path'] == path_select].iloc[0]
    
    # Create filter using boolean indexing
    mask = pd.Series(True, index=readme.index)
    filter_cols = ['Class', 'Experiment.Method', 'Group', 'Quant.Method', 
                  'logTransform', 'Normalized']
    
    for col in filter_cols:
        mask &= (readme[col] == params[col])
    
    return readme[mask]

def get_matching_path_with_pathology(readme, base_path, target_pathology):
    # Get parameters from base path
    params = readme[readme['Path'] == base_path].iloc[0]
    
    # Create filter using boolean indexing
    mask = pd.Series(True, index=readme.index)
    # Filter all columns except Path and Pathology
    filter_cols = ['Class', 'Experiment.Method', 'Group', 'Quant.Method', 
                  'logTransform', 'Normalized']
    
    # Match all parameters except Pathology
    for col in filter_cols:
        mask &= (readme[col] == params[col])
    
    # Add Pathology filter
    mask &= (readme['Pathology'] == target_pathology)
    
    # Get the matching path
    matching_rows = readme[mask]
    if matching_rows.shape[0] > 0:
        return matching_rows.iloc[0]['Path']
    return ""

def calc_diff(path1, path2,out_dir, gene_map=None):
    omics1=pd.read_csv(path1,sep="\t",header=0,index_col=0)
    omics2=pd.read_csv(path2,sep="\t",header=0,index_col=0)




    