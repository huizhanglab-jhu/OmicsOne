import matplotlib.pyplot as plt


import matplotlib.pyplot as plt


def plot_uniq(df,title,feature_type="Protein_Gene",figsize = (20,8), output_file=None):


    # Assuming the DataFrame has been correctly read into a variable called df


    # Plotting
    fig, ax = plt.subplots(figsize=figsize)
    # ax = plt.gca()
    

    # Bottom bar is shares
    ax.bar(df['SampleSet'], df['not_uniq'], label='Not unique')

    # Top bar is uniq, starting at the height of the shares
    ax.bar(df['SampleSet'], df['uniq'], bottom=df['not_uniq'], label='Unique')
    
    ax.bar(df['SampleSet'], df['common_100'], label='100%') 
    ax.bar(df['SampleSet'], df['common_50'], bottom=df['common_100'], label='>50%')

    # Rotation of x-tick labels
    # plt.xticks(rotation=90)

    # Adding legend
    ax.legend()

    # Adding labels and title
    ax.set_xlabel('SampleSet')
    ax.set_ylabel(f'{feature_type} (#)')
    ax.set_title(title)
    ax.yaxis.grid(True)
    ax.xaxis.grid(True, linestyle = '--')
    ax.set_axisbelow(True)
    
    max_num = 0
    for index,row in df.iterrows():
        num = row['uniq'] + row['not_uniq']
        if num > max_num:
            max_num = num
            
    plt.ylim(0,max_num*1.1)
    plt.xticks(rotation=45)
    # Display the plot
    plt.tight_layout()
    if output_file is None:
        return fig
    else:
        plt.savefig(output_file,bbox_inches='tight')

