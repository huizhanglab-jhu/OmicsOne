import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gseapy as gp


def omicsone_enrichr(genes,total_genes,outdir,**kwargs):
    gene_sets = kwargs.get('gene_sets',['KEGG_2021_Human'] )
    organism = kwargs.get('organism','human')
    fdr = kwargs.get('fdr',0.05)
    enr = gp.enrichr(gene_list=list(genes),
                                 background=len(total_genes),
                                 gene_sets=gene_sets,
                                 organism=organism,
                                 outdir = outdir
                                 )

    df = enr.res2d
    df = df[df['Adjusted P-value']<fdr]
    df2 = df.copy(deep=True)
    df2['-Log10(Adj.P)'] = df2['Adjusted P-value'].apply(lambda x: -1 * np.log10(x))
    df2 = df2.sort_values('-Log10(Adj.P)', ascending=False)
    return df2


def plot_enrichr_both(up_df,down_df,title,min_x=-10,max_x=10):
    rows = []
    for index,row in down_df.iterrows():
        d = row['Term']
        er = -1 * float(row['-Log10(Adj.P)'])
        c = 'Down'
        rows.append([d + " ",er,c])
    for index,row in up_df.iterrows():
        d = row['Term']
        er = float(row['-Log10(Adj.P)'])
        c = "Up"
        rows.append([" " + d,er,c])
    enrich_df = pd.DataFrame(rows,columns=['Term','-Log10(FDR)','Class']).sort_values('-Log10(FDR)',ascending=False)
    # enrich_df

    max_e = enrich_df['-Log10(FDR)'].max() + 1
    min_e = enrich_df['-Log10(FDR)'].min() - 1
    
    if max_e > max_x:
        max_x = max_e
        
    if min_e < min_x:
        min_x = min_e

    colors = ['red','blue']
    palette = sns.color_palette(colors)
    sns.set_style('white')
    fig, ax = plt.subplots(figsize=(4,6),dpi=100)
    ax1 = sns.barplot(data=enrich_df,x='-Log10(FDR)',y='Term',hue='Class', palette=palette,)
    
    n = 0
    for index,row in enrich_df.iterrows():
        if row['-Log10(FDR)'] > 0:
            ax.text(-0.5,n+0.5, row['Term'],horizontalalignment = 'right', verticalalignment = 'bottom')
        else:
            ax.text(0.5,n+0.5, row['Term'],horizontalalignment = 'left', verticalalignment = 'bottom')
        n += 1
    plt.yticks([])
    plt.legend(bbox_to_anchor=(1,0.85))
    # t = ax.set_xticklabels(['','3','2','1','0','1','2',''])
    plt.xlim(min_x,max_x)
    plt.title(title)
    plt.ylabel('')

    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    # plt.plot([0,0],[-3,15],color='k',linewidth=0.5)
    plt.axvline(x=0,color='k',lw=0.5)
    plt.axhline(y=up_df.shape[0] + down_df.shape[0],color='k',lw=0.5)
    # plt.plot([min_x,max_x],[15,15],color='k',linewidth=0.5)
    # plt.tick_params(top='off', left='off', right='off', labelleft='off', labelbottom='on')

    return fig