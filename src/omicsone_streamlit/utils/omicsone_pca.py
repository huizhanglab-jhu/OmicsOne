import streamlit as st
import numpy as np
import seaborn as sns
import os,sys,re
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
from sklearn import preprocessing

def plot_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radii.
    kwargs : dict
        Forwarded to `~matplotlib.patches.Ellipse`
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    # cov = np.cov(x, y)
    # pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # ell_radius_x = np.sqrt(1 + pearson)
    # ell_radius_y = np.sqrt(1 - pearson)
    # ellipse = Ellipse((np.mean(x), np.mean(y)),
    #                 width=ell_radius_x * 2 * n_std,
    #                 height=ell_radius_y * 2 * n_std,
    #                 angle=np.rad2deg(np.arctan2(*np.linalg.eig(cov)[1][:, 0][::-1])),
    #                 **kwargs)
    cov = np.cov(x, y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    ellipse = Ellipse(xy=(mean_x, mean_y),
                width=lambda_[0]*n_std*2, height=lambda_[1]*n_std*2,
                angle=np.rad2deg(np.arctan2(*v[:, 0][::-1])),
                    **kwargs)


    ellipse.set_edgecolor(kwargs.get('edgecolor', 'black'))
    ellipse.set_facecolor('none')
    ax.add_patch(ellipse)
    return ellipse



def calculate_pca_tables(data_df, component_path, score_path, corr_path,
                         variance_path, contribution_path, n_components=2):
    x = data_df.T.dropna()
    x = x.T
    features = x.columns.values
    samples = x.index
    data_scaled = pd.DataFrame(preprocessing.scale(x), columns=features, index=samples)
    pca = PCA(n_components=n_components)
    X = pca.fit_transform(data_scaled)
    component_names = ['PC-{}'.format(i + 1) for i in range(pca.components_.shape[0])]
    component_df = pd.DataFrame(pca.components_, columns=features,
                                index=component_names)
    contributions = []
    for index, row in component_df.iterrows():
        x = np.abs(np.array(list(row)))
        s = np.sum(x)
        new_row = x / s
        contributions.append(new_row)
    scaled_contributions = np.sum(pca.explained_variance_ratio_[:, None] * np.abs(
        np.array(contributions)
    ), axis=0)
    contributions.append(scaled_contributions)
    contribution_df = pd.DataFrame(contributions, columns=features, index=component_names + ['Scaled'])

    score_df = pd.DataFrame(X, index=samples, columns=component_names)
    corr_arr = pca.components_.T * np.sqrt(pca.explained_variance_)
    corr_df = pd.DataFrame(corr_arr, index=features, columns=component_names)
    component_df.to_csv(component_path, sep="\t", index=True)
    score_df.to_csv(score_path, sep="\t", index=True)
    corr_df.to_csv(corr_path, sep="\t", index=True)
    variance = pca.explained_variance_
    variance_ratio = pca.explained_variance_ratio_
    rows = []
    for i, j in zip(variance, variance_ratio):
        rows.append([i, j])
    variance_df = pd.DataFrame(rows, columns=['Variance', 'Variance_ratio'],
                               index=component_names)
    variance_df.to_csv(variance_path, sep="\t", index=True)
    contribution_df.to_csv(contribution_path, sep="\t", index=True)



def plot_pca(score_df, job_dir, figure_settings):

    palette = figure_settings['palette'] if 'palette' in figure_settings else 'Dark2'
    figure_width = figure_settings['width'] if 'width' in figure_settings else 6
    figure_height = figure_settings['height'] if 'height' in figure_settings else 6

    
    score_df = score_df.sort_values(by='Group')

    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(figure_width, figure_height))

    pc_x = 'PC-1'
    pc_y = 'PC-2'
    sns.scatterplot(data=score_df, x=pc_x, y=pc_y, hue='Group', palette=palette, 
                    ax=ax, s=50, alpha=0.7, legend='full')

    # 绘制椭圆
    n_groups = len(score_df['Group'].unique())
    for index, group in enumerate(sorted(list(score_df['Group'].unique()))):
        group_data = score_df[score_df['Group'] == group]
        if group_data.shape[0] == 1:
            continue
        try:
            plot_ellipse(group_data[pc_x], group_data[pc_y], ax=ax, 
                            edgecolor=sns.color_palette(palette, n_groups)[index])
        except Exception as e:
            st.error(f"Error plotting ellipse: {e}, {group}")

    ax.set_xlabel(pc_x)
    ax.set_ylabel(pc_y)

    plt.tight_layout()
    # plt.show()
    
    return fig



def pca_view(data_df, meta_df, meta_select, job_dir):


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
        selected_palette = st.selectbox("Select Color Palette", list(palettes.keys()) + ["Custom"])

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
                step=1  # Step size for increment/decrement
            )
        with size_cols[1]:
            figure_width= st.number_input(
                "Enter Figure size (width):",
                min_value=0,  # Minimum value (optional)
                max_value=100,  # Maximum value (optional)
                value=10,  # Default value
                step=1  # Step size for increment/decrement
            )

        figure_settings = {
            'width': figure_width,
            'height': figure_height,
            'palette': custom_palette,
        }

    with left_col:
        with st.spinner('Running... Please wait...'):
            component_path = os.path.join(job_dir, 'pca_component.tsv')
            score_path = os.path.join(job_dir, 'pca_score.tsv')
            corr_path = os.path.join(job_dir, 'pca_corr.tsv')
            variance_path = os.path.join(job_dir, 'pca_variance.tsv')
            contribution_path = os.path.join(job_dir, 'pca_contribution.tsv')

            
            data = data_df.iloc[:, :-1]
            calculate_pca_tables(data, component_path, score_path, corr_path,
                         variance_path, contribution_path, n_components=2)
            score_df = pd.read_csv(score_path, sep='\t', index_col=0)
            score_df['Group'] = data_df['Group']
            fig = plot_pca(score_df, job_dir, figure_settings)
            st.pyplot(fig)
