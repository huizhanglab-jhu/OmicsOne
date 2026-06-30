import statsmodels.stats.multitest as smt
from collections import defaultdict

import pandas as pd
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from scipy.cluster.hierarchy import linkage
import seaborn as sns
from matplotlib.colors import rgb2hex, colorConverter
import os,sys,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
from matplotlib.pyplot import figure, show
from scipy.stats import spearmanr, variation,zscore
import io
# import cgi
import base64
from collections import OrderedDict


def labels2cmap(labels, sample_info, color_sys='hls'):
    import seaborn as sns
    si = sample_info.replace(np.nan, 'N/A')
    cmap = {}
    for label in labels:
        x = si.loc[:, label].unique()
        # print(x)
        cmap[label] = dict(zip(x, sns.color_palette(color_sys, len(list(x)))))
    print(cmap)
    return cmap


# def labels2colors(labels,sample_info,key_col='SPL',color_sys='hls'):
#     cmap = labels2cmap(labels,sample_info,color_sys)
#     si = sample_info.loc[:, [key_col] + labels]

def labels2colors(labels, sample_info, key_col='SPL', color_sys='hls'):
    cmap = labels2cmap(labels, sample_info, color_sys)
    si = sample_info.loc[:, [key_col] + labels]
    si = si.set_index(key_col)
    for label in labels:
        si[label] = si[label].map(cmap[label])
    return si
    



def newclustermap2(df,name,out_dir):
    df=pd.DataFrame(zscore(df,axis=1,ddof=1),columns=df.columns, index=df.index)
    # link = linkage(df.T, metric='correlation', method='ward')
    # Calculate correlation distance matrix and convert to Euclidean
    corr_matrix = np.corrcoef(df.T)
    # Convert correlation to distance (1 - correlation)
    dist_matrix = 1 - corr_matrix
    # Convert to Euclidean distance
    dist_matrix_euclidean = np.sqrt(2 * dist_matrix)
    
    # Use the Euclidean distances with Ward's method
    link = linkage(dist_matrix_euclidean, method='ward')
    
    den = dendrogram(link, labels=df.columns)
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    
    corr=[]
    # cluster_classes = Clusters()
    for c, l in cluster_idxs.items():
        i_l = [den['ivl'][i] for i in l]
        #cluster_classes[c] = i_l
        for x in range(len(i_l)): 
            corr.append([c,i_l[x]])
            
    corr2=pd.DataFrame(corr) 
    corr2.columns=["group","SPL"]
    flatui = ["g", "r", "c", "m", "y", "k","b"]
    sns.palplot(sns.color_palette(flatui))
    mycolors2 = labels2cmap(["group"], corr2, color_sys=flatui)
    col_colors2a = labels2colors(["group"],corr2,color_sys=flatui)
    col_colors2a =  col_colors2a.loc[~ col_colors2a.index.duplicated(keep='first')]
    g = sns.clustermap(df, z_score=0,vmax=2,vmin=-2,cmap='RdBu_r',col_colors=col_colors2a, 
                   col_linkage=link,figsize=(5, 5),xticklabels=False, yticklabels=False
                  )
    # Add axis labels
    g.ax_heatmap.set_xlabel('Samples')
    g.ax_heatmap.set_ylabel('Gene')
    fig = g.fig
    corr3=corr2
    corr3.set_index('SPL', inplace=True,drop=True)
    corr3=corr3.loc[~ corr3.index.duplicated(keep='first')]
    plt.savefig(os.path.join(out_dir,f"{name}.png"),dpi=100,bbox_inches='tight')
    return corr3, fig


def calculate_gene_high_and_low(omics1a, omics2a, out_dir):
    
    corr=[]
    for index, row in omics1a.iterrows():
        corr1,p=spearmanr(omics1a.loc[index,:].T,omics2a.loc[index,:].T)
        cv1 = variation(omics1a.loc[index,:], axis=0)
        cv2 = variation(omics1a.loc[index,:], axis=0)
        corr.append([index,corr1,p,cv1,cv2])
    #     print([index,corr1,p])
    df_corrtumorRNAproteinGene=pd.DataFrame(np.asarray(corr))
    df_corrtumorRNAproteinGene.apply(pd.to_numeric, errors='ignore', downcast='float')
    df_corrtumorRNAproteinGene.columns=["Gene","Tumor Spearman Correlation","Tumor P","CV1","CV2"]
    df_corrtumorRNAproteinGene["Tumor Spearman Correlation"]=df_corrtumorRNAproteinGene["Tumor Spearman Correlation"].astype(float)
    df_corrtumorRNAproteinGene.set_index('Gene', inplace=True)
    df_corrtumorRNAproteinGene['Tumor BH adjusted P'] =  smt.multipletests(df_corrtumorRNAproteinGene['Tumor P'].astype(float), method='fdr_bh')[1]
    # # mRNA-protein correlation
    if(len(df_corrtumorRNAproteinGene.index)>1000):
        cvhigh=df_corrtumorRNAproteinGene.astype(float).nlargest(1000,'CV1').index
    else:
        cvhigh=df_corrtumorRNAproteinGene.astype(float).nlargest(round(len(df_corrtumorRNAproteinGene.index)),'CV1').index
    genehigh=df_corrtumorRNAproteinGene.astype(float).nlargest(round(len(df_corrtumorRNAproteinGene.index)/20),'Tumor Spearman Correlation').index
    genelow=df_corrtumorRNAproteinGene.astype(float).nsmallest(round(len(df_corrtumorRNAproteinGene.index)/20),'Tumor Spearman Correlation').index
    genehigh=set(cvhigh).intersection(set(genehigh))
    genelow=set(genelow).intersection(set(cvhigh))

    genehigh = sorted(list(genehigh))
    genelow = sorted(list(genelow))

    return genehigh, genelow