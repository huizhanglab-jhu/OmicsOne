

import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import os,sys,re
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_tsne_table(data_df, perplexity=30, n_components=2, random_state=42, 
                         output_path="tsne_results.tsv", corr_path="tsne_corr.tsv"):
    """
    Perform t-SNE dimensionality reduction and save results to a table.

    Parameters:
        data_df (pd.DataFrame): Input data (features as columns, samples as rows).
        perplexity (int): Perplexity parameter for t-SNE.
        n_components (int): Number of dimensions for t-SNE embedding.
        random_state (int): Random state for reproducibility.
        output_path (str): Path to save the t-SNE result table.

    Returns:
        tsne_df (pd.DataFrame): t-SNE results as a DataFrame.
    """
    # Preprocess the data: Standardize features
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_df)
    
    # Perform t-SNE
    tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=random_state)
    tsne_results = tsne.fit_transform(data_scaled)
    
    # Create a DataFrame for the t-SNE results
    tsne_columns = [f"t-SNE-{i+1}" for i in range(n_components)]
    tsne_df = pd.DataFrame(tsne_results, index=data_df.index, columns=tsne_columns)
    
    # Save to file
    tsne_df.to_csv(output_path, sep="\t")
    print(f"t-SNE results saved to: {output_path}")

    # Calculate correlations between original features and t-SNE components
    corr_arr = np.corrcoef(data_scaled.T, tsne_results.T)[:data_scaled.shape[1], data_scaled.shape[1]:]
    corr_df = pd.DataFrame(corr_arr, index=data_df.columns, columns=tsne_columns)
    
    # Save t-SNE results and correlations to files
    # tsne_df.to_csv(output_path, sep="\t")
    corr_df.to_csv(corr_path, sep="\t")
    
    return tsne_df, corr_df

def plot_tsne(tsne_df, labels=None, figsize=(10, 8), cmap="viridis"):
    """
    Plot the t-SNE embedding.

    Parameters:
        tsne_df (pd.DataFrame): t-SNE results (e.g., from calculate_tsne_table).
        labels (pd.Series or list): Optional labels for coloring points.
        figsize (tuple): Figure size for the plot.
        cmap (str): Colormap for coloring points.
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if labels is not None:
        sns.scatterplot(
            x=tsne_df.iloc[:, 0],
            y=tsne_df.iloc[:, 1],
            hue=labels,
            palette=cmap,
            legend="full",
            s=50,
            ax=ax
        )
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        sns.scatterplot(
            x=tsne_df.iloc[:, 0],
            y=tsne_df.iloc[:, 1],
            legend=None,
            s=50,
            ax=ax
        )
    
    ax.set_xlabel(tsne_df.columns[0])
    ax.set_ylabel(tsne_df.columns[1])
    ax.set_title("t-SNE Embedding")
    fig.tight_layout()
    # plt.show()
    return fig,ax



def tsne_view(data_df, meta_df, meta_select, job_dir):


    left_col, right_col = st.columns(2)

    meta_samples = meta_df.index.tolist()
    data_samples = list(data_df.columns.values)

    samples = [i for i in data_samples if i in meta_samples]
    sample_map = dict(zip(meta_samples, list(meta_df[meta_select])))
    
    # only keep samples in both data and meta
    data_df = data_df[samples]

    # preprocess data
    # convert None to nan, and convert values in dataframe to float
    data_df = data_df.replace('None',np.nan) \
        .applymap(float) \
        .dropna()

    # Check for NaNs
    if data_df.isnull().values.any():
        st.error("Data contains NaNs. Please clean your data.")
        return

    # Check for infinite values
    if np.isinf(data_df).values.any():
        st.error("Data contains infinite values. Please clean your data.")
        return
    
    data_df = data_df.T
    data_df['Group'] = [sample_map[index] for index,row in data_df.iterrows()]
        
    with right_col:
        # Define color palettes
        uniq_groups = data_df['Group'].unique()
        # Define predefined color palettes
        palettes = {
            'Dark2': sns.color_palette('Dark2', n_colors=len(uniq_groups)),
            'Set1': sns.color_palette('Set1', n_colors=len(uniq_groups)),
            'Pastel1': sns.color_palette('Pastel1', n_colors=len(uniq_groups)),
            'Paired': sns.color_palette('Paired', n_colors=len(uniq_groups)),
        }

        # Streamlit color palette selector
        selected_palette = st.selectbox("Select Color Palette", list(palettes.keys()) + ["Custom"], key='tsne_palette')

        # Initialize the color palette
        if selected_palette == "Custom":
            # Allow user to input a custom palette
            custom_colors_input = st.text_input("Enter custom colors (comma-separated)", 
                                                "#1f77b4, #d62728")
            # Convert input string to a list of colors
            custom_palette = [color.strip() for color in custom_colors_input.split(',')]
        else:
            # Use the selected predefined palette
            custom_palette = palettes[selected_palette]

        size_cols = st.columns([1,1])
        with size_cols[0]:
            figure_height= st.number_input(
                "Enter Figure size (height):",
                min_value=0,  # Minimum value (optional)
                max_value=100,  # Maximum value (optional)
                value=10,  # Default value
                step=1,  # Step size for increment/decrement
                key='tsne_figure_height'
            )
        with size_cols[1]:
            figure_width= st.number_input(
                "Enter Figure size (width):",
                min_value=0,  # Minimum value (optional)
                max_value=100,  # Maximum value (optional)
                value=10,  # Default value
                step=1 , # Step size for increment/decrement
                key='tsne_figure_width'
            )

        figure_settings = {
            'width': figure_width,
            'height': figure_height,
            'palette': custom_palette,
        }

    with left_col:
        with st.spinner('Running... Please wait...'):
            tsne_score_path = os.path.join(job_dir, 'tsne_score.tsv')
            corr_path = os.path.join(job_dir, 'tsne_corr.tsv')
            
            data = data_df.iloc[:, :-1]

            tsne_df , corr_df = calculate_tsne_table(data, perplexity=30, n_components=2, 
                                                     output_path=tsne_score_path, corr_path=corr_path)
            labels = data_df['Group']
            fig,ax = plot_tsne(tsne_df, labels=labels, cmap=custom_palette, figsize=(figure_width, figure_height))
            # fig = plot_pca(score_df, job_dir, figure_settings)
            st.pyplot(fig)
