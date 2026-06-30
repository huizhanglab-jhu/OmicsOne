import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

def plot_volcano(diff_df,s=5, hue='Significance', **kwargs):

    x_min = np.min(diff_df['Log2FC(median)'])
    x_max = np.max(diff_df['Log2FC(median)'])
    abs_x_max = max(abs(x_min), abs(x_max))
    x_min = -abs_x_max
    x_max = abs_x_max
    
    y_max = np.max(diff_df['-Log10(FDR)'])

    color_dict = {
        'NS': 'gray',
        'D': 'skyblue',
        'U': 'salmon',
        'S-D': 'blue',
        'S-U': 'red'
    }

    if hue == "GlycanType":
        color_dict = {
            "HM": "#1f77b4",      # Blue
            "only_F": "#ff7f0e",  # Orange
            "only_S": "#2ca02c",  # Green
            "F+S": "#d62728",     # Red
            "Other": "#9467bd"    # Purple
        }
    

    fig, ax = plt.subplots(figsize=(8,8))
    scatter = sns.scatterplot(data = diff_df, x = 'Log2FC(median)', y = '-Log10(FDR)', 
                    hue = hue,s=s, palette=color_dict,edgecolor='none')
    
    if 'annotations' in kwargs:


        pos_map = dict()
        for index,row in diff_df.iterrows():
            if index in kwargs['annotations'] :       
                x = float(row["Log2FC(median)"])
                y = float(row["-Log10(FDR)"])
                pos_map[index] = (x,y)
                # print(key,glycan,index[3],index[2])

        y_values = [pos_map[key][1] for key in pos_map]
        y_values_sorted = sorted(y_values)
        # print(y_values_sorted)
        # print(np.min(y_values_sorted),np.max(y_values_sorted))

        old_y_max = np.max(y_values_sorted)
        old_y_min = np.min(y_values_sorted)
        # print(old_y_max,old_y_min)

        new_y_max = old_y_max + 1
        new_y_min = old_y_min - 5
        new_y_values = np.linspace(new_y_min, new_y_max, len(y_values_sorted))
        # print(new_y_values)
        y_mapping = {old_y: new_y for old_y, new_y in zip(y_values_sorted, new_y_values)}
        # print(y_mapping)

        old_y_range = old_y_max - old_y_min
        new_y_range = new_y_max - new_y_min

        new_pos_map = dict()

        for key in pos_map:
            x,y = pos_map[key]
            new_y = y_mapping[y]
            new_pos_map[key] = (x-2 if x < 0 else x + 2, new_y)


        for key in new_pos_map:
                x_text,y_text = new_pos_map[key]
                x,y = pos_map[key]
                data_type =  kwargs.get('data_type', 'glycosite') 
                if data_type == "glycosite":
                    protein, gene, site, glycan = key.split('_')
                    glycan = re.sub("[FSG]0","",glycan)
                    new_key = f"{gene}_{site}_{glycan}"
                elif data_type == "phosphosite":
                    gene_id, protein_id, site, sequence, _, gene = key.split('|')
                    new_key = f"{gene}_{site}"
                else:
                    new_key = key
                
                if x < 0:
                    ax.text(x_text, y_text, new_key, ha='right',
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", 
                                facecolor="white", alpha=0.5))
                else:
                    ax.text(x_text, y_text, new_key, ha='left',
                            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", 
                                    facecolor="white", alpha=0.5))
                ax.annotate("",xy=(x,y),xytext=(x_text, y_text),ha='center',
                    bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", 
                              facecolor="white", alpha=0.5),
                              arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2", color="black", alpha=0.2))


    
    # Add horizontal line at -log10(0.01)
    plt.axhline(y=-np.log10(0.01), color='black', linestyle='--', alpha=0.5)

    # Add vertical lines at log2(2) and -log2(2)
    log2fc_threshold = kwargs.get('log2fc_threshold', np.log2(2))
    plt.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-1 * log2fc_threshold, color='black', linestyle='--', alpha=0.5)
   
    # Move the legend to the bottom right
    # ax.legend( loc='lower right')
    # Move the legend to the bottom right and set marker size and title
    from matplotlib.legend_handler import HandlerPathCollection
    legend = ax.legend(loc='lower right', title='Significance', prop={'size': 10},markerscale=3)
    # for handle in legend.legendHandles:
    #     handle._sizes = [30]  # Set marker size

    
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)

    if xlim is not None:
        plt.xlim(xlim)
    else:
        plt.xlim(x_min-11-0.1,x_max+11+0.1)

    if ylim is not None:
        plt.ylim(0,y_max + 0.1)
    else:
        plt.ylim(ylim)


    return fig
