import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
    
import math

import pandas as pd

from matplotlib.patches import Rectangle

from utils.genome import cis_trans_corr_count
import numpy as np

def scatter_wrapper(ax, df, x='cnv', y='protein', s=None, c=None, marker=None, cmap=None, norm=None, vmin=None,
                    vmax=None, alpha=None, linewidths=None,  edgecolors=None,
                    **kwargs):
    # print(df.head(2))
    ax.scatter(np.array(df[x], dtype=float),
               np.array(df[y], dtype=float),
               s=s, c=c, marker=marker, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, alpha=alpha, linewidths=linewidths,
                edgecolors=edgecolors)
    


# Function to draw a rectangle
def draw_rectangle(ax, startx, starty, width, height, color='blue', label=None):
    rect = Rectangle((startx, starty), width, height, color=color, label=label)
    ax.add_patch(rect)

def draw_corr_count(ax,final_corr_count,chr_bounds, total_bp , top_n = 5,):

    top_arms = set(final_corr_count.sort_values('count.sum',ascending=False).head(top_n)['cnv.arm'])

    # Iterate over rows to plot rectangles
    for _, row in final_corr_count.iterrows():
        
        # Base rectangle (pos.cis)
        startx = row['cnv.arm.start']
        starty = 0
        height = row['pos.cis.count']
        width = row['cnv.arm.end'] - row['cnv.arm.start']
        if height > 0:
            draw_rectangle(ax, startx, starty, width,height, color='pink')

        # stacked (pos.trans)
        starty = row['pos.cis.count']
        height = row['pos.trans.count']
        if height > 0:
            draw_rectangle(ax, startx, starty, width,height,color='red')

        starty = row['neg.cis.count'] * (-1)
        height = row['neg.cis.count'] 
        if height > 0:
            draw_rectangle(ax, startx, starty, width,height, color='skyblue')

        starty = (row['neg.cis.count']   + row['neg.trans.count']) * (-1)
        height = row['neg.trans.count']
        if height > 0:
            draw_rectangle(ax, startx, starty, width,height, color='blue')

        height = row['pos.cis.count'] + row['pos.trans.count']

        if row['cnv.arm'] in top_arms:
            ax.text(startx, height, row['cnv.arm'])
    
    # print("DEBUG: final_corr_count")
    # print(final_corr_count['neg.cis.count'], final_corr_count['neg.trans.count'])
    # print(final_corr_count['pos.cis.count'], final_corr_count['pos.trans.count'])
    ax.set_ylim(np.max([final_corr_count['neg.cis.count'], final_corr_count['neg.trans.count']]) * (-1) * 1.1 - 1000,
                np.max([final_corr_count['pos.cis.count'], final_corr_count['pos.trans.count']]) * 1.5)
    ax.set_xlim(0,total_bp)
    # ax.set_xticklabels([])
    ax.set_xticks(chr_bounds)
    ax.set_xticklabels([])
    ax.tick_params(direction='out')


def draw_corr(corr_df, **kwargs):
    
    fig = plt.figure(figsize=(5, 10),dpi=100)
    # add_axes([left_x,left_y,width,height])
    ax1 = fig.add_axes([0.05, 0.1, 0.9, 0.02])
    ax2 = fig.add_axes([0.05, 0.15, 0.9, 0.1])
    ax3 = fig.add_axes([0.05, 0.4, 0.9, 0.4])
    ax4 = fig.add_axes([0.05, 0.28, 0.9, 0.1])
    
    corr_with = kwargs.get('corr_with',None)
    if corr_with is None:
        print("corr_with is None")
        return None
    
    chr_bounds = kwargs.get('chr_bounds',None)
    if chr_bounds is None:
        print("chr_bounds is None")
        return None
    
    chr_names = kwargs.get('chr_names',None)
    if chr_names is None:
        print("chr_names is None")
        return None
    
    chr_ticks = kwargs.get('chr_ticks',None)
    if chr_ticks is None:
        print("chr_ticks is None")
        return None
    
    chr_names2 = kwargs.get('chr_names2',None)
    if chr_names2 is None:
        print("chr_names2 is None")
        return None
    
    genes = kwargs.get('genes',None)
    if genes is None:
        print("genes is None")
        return None
    
    samples = kwargs.get('samples',None)
    if samples is None:
        print("samples is None")
        return None
    
    gistic_data = kwargs.get('gistic_data',None)
    if gistic_data is None:
        print("gistic_data is None")
        return None
    
    total_bp = kwargs.get('total_bp',None)
    if total_bp is None:
        print("total_bp is None")
        return None
    
    gene_map = kwargs.get('gene_map',None)
    if gene_map is None:
        print("gene_map is None")
        return None
    
    start_map = kwargs.get('start_map',None)
    if start_map is None:
        print("start_map is None")
        return None
    
    cytoband_d = kwargs.get('cytoband_d',None)
    if cytoband_d is None:
        print("cytoband_d is None")
        return None 

    # plot the bottom chromosome
    color_list = ['black' if i % 2 ==0 else 'white' for i in range(24)]
    cmap = mpl.colors.ListedColormap(color_list)
    norm = mpl.colors.BoundaryNorm(chr_bounds, cmap.N)
    cb2 = mpl.colorbar.ColorbarBase(ax1,
                                    cmap=cmap,
                                    norm=norm,
                                    ticks=chr_ticks,
                                    spacing='proportional',
                                    orientation='horizontal')
    ax1.set_xticklabels(chr_names,  rotation=0 )
    plt.setp( ax1.xaxis.get_majorticklabels(), ha="center" )

    if corr_with == 'RNA':
        scatter_wrapper(ax3,corr_df,'cnv.location','rna.location',
                            c = corr_df['Correlation'], cmap='bwr', marker="_",s =10 ,alpha=0.05)
    elif corr_with == 'Protein':
        scatter_wrapper(ax3,corr_df,'cnv.location','protein.location',
                            c = corr_df['Correlation'], cmap='bwr', marker="_",s =10 ,alpha=0.05)
    ax3.set_xlim([0, 3095677412])
    ax3.set_ylim([0, 3095677412])
    ax3.set_xticks(chr_bounds)
    ax3.set_xticklabels([])
    ax3.set_yticks(chr_bounds)
    ax3.set_yticklabels(['']+chr_names2)
    
    ax2.tick_params(direction='out')

    ax3.patch.set_visible(False)
    ax3.grid(alpha=0.2)

    final_corr_count = cis_trans_corr_count(corr_df, start_map, cytoband_d, corr_with)
    # print(final_corr_count.head(2))
    draw_corr_count(ax2,final_corr_count,chr_bounds,total_bp, top_n = 8)
    ax2.grid(alpha=0.2)
    
    from utils.genome import calc_gistic
    
    

    line_segments_amp, line_segments_del = calc_gistic(gistic_data,genes, samples, gene_map, start_map)
    
    # print(line_segments_amp)
    # print(line_segments_del)

    amp_lc = LineCollection(line_segments_amp, colors="red", linewidths=2)
    del_lc = LineCollection(line_segments_del, colors="blue", linewidths=2)
    ax4.add_collection(amp_lc)
    ax4.add_collection(del_lc)
    ax4.set_xlim(0,total_bp)
    ax4.set_ylim(-100,100)
    ax4.set_xticks(chr_bounds)
    ax4.set_xticklabels([])
    ax4.grid(alpha=0.2)
    # ax4.add_collection(lc2)
    # ax4.autoscale()
    ax4.autoscale_view()
    
    return fig
