

import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import os,sys,re

def calculate_umap(data_df, n_neighbors=15, n_components=2, random_state=42, 
                   embedding_path='umap_embedding.tsv', corr_path='umap_corr.tsv'):
    """
    Perform UMAP dimensionality reduction, calculate correlations, and save results.

    Parameters:
        data_df (pd.DataFrame): Input data (features as columns, samples as rows).
        n_neighbors (int): Number of neighbors for UMAP.
        n_components (int): Number of dimensions for UMAP embedding.
        random_state (int): Random state for reproducibility.
        output_prefix (str): Prefix for saving result files.

    Returns:
        embedding_df (pd.DataFrame): UMAP embedding results.
        corr_df (pd.DataFrame): Correlations between original features and UMAP components.
    """

    # Data preprocessing: Scale features
    samples = data_df.index
    features = data_df.columns
    data_scaled = pd.DataFrame(preprocessing.scale(data_df), columns=features, index=samples)
    
    # Import lazily so Streamlit startup and non-UMAP pages do not pay the numba import cost.
    from umap.umap_ import UMAP
    umap_reducer = UMAP(n_neighbors=n_neighbors, n_components=n_components, random_state=random_state)
    embedding = umap_reducer.fit_transform(data_scaled)
    
    # Save UMAP embedding
    component_names = [f"UMAP-{i + 1}" for i in range(n_components)]
    embedding_df = pd.DataFrame(embedding, index=samples, columns=component_names)
    embedding_df.to_csv(embedding_path, sep="\t")
    print(f"UMAP embedding saved to: {embedding_path}")
    
    # Calculate correlations
    corr_arr = np.corrcoef(data_scaled.T, embedding.T)[:len(features), len(features):]
    corr_df = pd.DataFrame(corr_arr, index=features, columns=component_names)
    corr_df.to_csv(corr_path, sep="\t")
    print(f"Feature correlations saved to: {corr_path}")
    
    return embedding_df, corr_df

def plot_umap(embedding_df, labels=None, figsize=(10, 8), cmap="viridis"):
    """
    Plot the UMAP embedding using fig and ax.

    Parameters:
        embedding_df (pd.DataFrame): UMAP embedding results (e.g., from calculate_umap).
        labels (pd.Series or list): Optional labels for coloring points.
        figsize (tuple): Figure size for the plot.
        cmap (str): Colormap for coloring points.
    """
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    if labels is not None:
        sns.scatterplot(
            x=embedding_df.iloc[:, 0],
            y=embedding_df.iloc[:, 1],
            hue=labels,
            palette=cmap,
            legend="full",
            s=50,
            ax=ax  # Pass the axis object
        )
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        sns.scatterplot(
            x=embedding_df.iloc[:, 0],
            y=embedding_df.iloc[:, 1],
            legend=None,
            s=50,
            ax=ax  # Pass the axis object
        )
    
    # Set labels and title
    ax.set_xlabel(embedding_df.columns[0])
    ax.set_ylabel(embedding_df.columns[1])
    ax.set_title("UMAP Embedding")
    
    # Adjust layout
    fig.tight_layout()
    
    return fig, ax





def umap_view(data_df, meta_df, meta_select, job_dir):


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
        selected_palette = st.selectbox("Select Color Palette", list(palettes.keys()) + ["Custom"], key='umap_palette')

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
                key='umap_figure_height'
            )
        with size_cols[1]:
            figure_width= st.number_input(
                "Enter Figure size (width):",
                min_value=0,  # Minimum value (optional)
                max_value=100,  # Maximum value (optional)
                value=10,  # Default value
                step=1 , # Step size for increment/decrement
                key='umap_figure_width'
            )

        figure_settings = {
            'width': figure_width,
            'height': figure_height,
            'palette': custom_palette,
        }

    with left_col:
        with st.spinner('Running... Please wait...'):
            embedding_path = os.path.join(job_dir, 'umap_embedding.tsv')
            corr_path = os.path.join(job_dir, 'umap_corr.tsv')
            
            data = data_df.iloc[:, :-1]
            embedding_df, corr_df = calculate_umap(data, n_neighbors=15, n_components=2, random_state=42, 
                           embedding_path=embedding_path, corr_path=corr_path)
            labels = data_df['Group']
            fig,ax = plot_umap(embedding_df, labels=labels, cmap=custom_palette, figsize=(figure_width, figure_height))
            # fig = plot_pca(score_df, job_dir, figure_settings)
            st.pyplot(fig)
