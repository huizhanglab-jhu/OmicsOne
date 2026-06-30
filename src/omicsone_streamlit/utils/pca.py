import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from .utils import print_result

np.random.seed(0)
from sklearn import preprocessing
import matplotlib.pyplot as plt
import plotly.express as px


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




@print_result
def pca_standalone(X, sample_info, param={}):
    group = param.get('group', 'Pathological_Status(Categorical)')
    # X.to_csv("exported_group_dfX.csv", index=True)
    old_group_name = list(X["group"])[0]
    group = list(X["group"])[0]
    for i in ['(Categorical)', "(Numeric)", "_"]:
        if i in group:
            if i == "_":
                group = group.replace(i, " ")
            else:
                group = group.replace(i, "")
    sample_info.columns = [group if i == old_group_name else i for i in sample_info.columns]
    X.drop(["group"], axis=1, inplace=True)
    imputation = param.get('imputation', None)

    pca_df = X.copy(deep=True)
    figure_height = param.get('figure_height', 500)
    if imputation == 'dropna':
        pca_df.dropna(inplace=True)
    colors = param.get('colors', ["red", "blue", "orange", "grey", "pink", "black",
                                  "navy", "purple", "green", "yellow"])
    fig = two_demensional_reduction(pca_df, sample_info, group=group,
                                    figure_height=figure_height,
                                    colors=colors, method='PCA')
    return fig


def two_demensional_reduction(orig_df, orig_si, group, method, colors, figure_height=None):
    si = orig_si[orig_si.index.isin(orig_df.columns.values)]
    si = si.applymap(lambda i: 'NA' if pd.isnull(i) else i)
    df = orig_df.copy(deep=True)
    df = df[[i for i in df.columns.values if i in si.index]]
    features = df.index
    X = df.T
    samples = X.index
    data_scaled = pd.DataFrame(preprocessing.scale(X), columns=features, index=samples)
    pca_var, pca_X = run_pca(data_scaled)
    si = si.loc[pca_X.index, [group]]
    pca_X[group] = si[group].values
    fig = draw_scatter_plot_new(pca_X, group, colors, figure_height)
    return fig


def run_pca(data_scaled):
    pca = PCA(n_components=2)
    X = pca.fit_transform(data_scaled)
    df_pca_var = pd.DataFrame(pca.components_, columns=data_scaled.columns,
                              index=['PC-1', 'PC-2'])
    pca_X = pd.DataFrame(X, index=data_scaled.index, columns=['PC-1', 'PC-2'])
    pca_var = df_pca_var.T
    return pca_var, pca_X


def plot_pca_scores(pca_score_df, pc_x='PC-1', pc_y='PC-2', height=500, width=500, group='Group', ellipse=True,
                    labels=[], title='',colors=None):
    groups = [str(i) for i in pca_score_df[group].unique()]
    colors = ["blue", "red", "orange", "grey", "pink", "black", "navy", "purple", "green", "yellow","dimgray",
              "rosybrown",'salmon','coral','saddlebrown','teal','cyan','violet','crimson'] if colors is None else colors
    color_discrete_map = {}
    for i, g in zip(range(len(groups)), sorted(groups)):
        color_discrete_map[g] = colors[i]
    pca_score_df['Sample'] = pca_score_df.index
    if labels is not None and len(labels) == 0:
        fig = px.scatter(pca_score_df, x=pc_x, y=pc_y, color=group, title=title,
                         color_discrete_map=color_discrete_map,

                         hover_data={
                             'Sample': True,
                         }
                         )
    else:
        fig = px.scatter(pca_score_df, x=pc_x, y=pc_y, color=group, title=title,
                         color_discrete_map=color_discrete_map,
                         text=pca_score_df.index if labels is None else labels,
                         hover_data={
                             'Sample': True,
                         }
                         )
    fig.for_each_trace(lambda t: t.update(textfont_color=t.marker.color, textposition='top center'))

    fig.update_traces(textposition='top center')
    axis_template = dict(
        linecolor='black',
        showline=True,
        showgrid=False,
        # mirror=True,
    )
    fig.update_layout(clickmode='event+select',
                      font_color='black',
                      legend=dict(
                          orientation="h",
                          yanchor="bottom",
                          y=1.02,
                          xanchor="right",
                          x=1,
                          # title=title
                      ),
                      xaxis=axis_template,
                      yaxis=axis_template,
                      margin=dict(l=20, r=20, t=20, b=20),
                      height=height,
                      width=width,
                      paper_bgcolor="White",
                      plot_bgcolor='White')

    if ellipse:
        for i in groups:
            # print(i)
            ps = pca_score_df[pca_score_df[group] == i][[pc_x, pc_y]]
            # print(ps)
            if ps.shape[0] > 2:
                fig.add_shape(type='path', path=plot_point_cov(ps),
                              line={'dash': 'dot'}, line_color="black")
    return fig


def plot_pca_contributions(pca_contribution_df, pc_name='Scaled', height=500):
    row = pca_contribution_df[pca_contribution_df.index == pc_name].iloc[0]
    values = np.array(row, dtype=float)
    orders = np.argsort(values)[::-1]
    size = len(values) * 0.1
    size = np.min([10, int(size)])
    orders = orders[:int(size)]
    values = [values[i] * 100 for i in orders]
    features = pca_contribution_df.columns.values
    names = [features[i] for i in orders]
    cols = [values, names]
    df = pd.DataFrame(cols).T
    df.columns = ['Contribution', 'Feature']
    df.index = names
    fig = px.bar(df, x='Feature', y='Contribution')
    fig.update_layout(clickmode='event+select',
                      margin=dict(l=20, r=20, t=20, b=20), height=height,
                      paper_bgcolor="LightSteelBlue")
    return fig


def plot_pca_variance(pca_variance_df, pca_contribution_df, height=500):
    rows = []
    for index, row in pca_variance_df.iterrows():
        new_row = [index, row['Variance_ratio'] * 100]
        rows.append(new_row)
    cols = pca_contribution_df.columns.values
    d = dict()
    for index, row in pca_contribution_df.iterrows():
        values = np.abs(np.array(list(row)))
        orders = np.argsort(values)[::-1][:10]
        features = [cols[i] for i in orders]
        d[index] = features
    df = pd.DataFrame(rows, columns=['PC', 'Variance_ratio'])

    rows = []
    for index, row in df.iterrows():
        pc = row['PC']
        for i in range(1, len(d[pc]) + 1):
            x = d[pc]
            row['Feature{}'.format(i)] = x[i - 1]
        rows.append(row)
    df = pd.DataFrame(rows)
    # print(df)
    hover_data = {}
    for i in range(1, 11):
        hover_data['Feature{}'.format(i)] = True
    fig = px.bar(df, x='PC', y='Variance_ratio',
                 hover_data=hover_data)
    axis_template = dict(
        linecolor='black',
        showline=True,
        # mirror=True,
    )
    fig.update_layout(clickmode='event+select',
                      xaxis=axis_template,
                      yaxis=axis_template,
                      xaxis_title='PCs',
                      yaxis_title='Variance ratio (%)',
                      margin=dict(l=20, r=20, t=20, b=20), height=height,
                      # mode="markers+lines",
                      plot_bgcolor="White",
                      paper_bgcolor="White")
    return fig


def draw_scatter_plot_new(X_df, group, colors, figure_height):
    fig = px.scatter(X_df, x="PC-1", y="PC-2", color=group, title="PCA")
    fig.update_layout(clickmode='event+select',
                      margin=dict(l=20, r=20, t=20, b=20), height=500,
                      paper_bgcolor="LightSteelBlue")

    groups = X_df[group].unique()
    # print(groups)
    for i in groups:
        ps = X_df[X_df[group] == i][["PC-1", "PC-2"]]
        # print("ps",ps)
        fig.add_shape(type='path', path=plot_point_cov(ps),
                      line={'dash': 'dot'}, line_color="black")

    return fig


def plot_point_cov(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).
    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.
    Returns
    -------
        A matplotlib ellipse artist
    """
    # pos = points.mean(axis=0)
    # print("pos::",pos)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse_1(cov, points, **kwargs)


def plot_cov_ellipse_1(cov, ps, n_std=1.96, size=500, **kwargs):
    import numpy as np

    """reference:https://gist.github.com/dpfoose/38ca2f5aee2aea175ecc6e599ca6e973
    https://community.plotly.com/t/how-to-draw-ellipse-on-top-of-scatter-plot/36576
    https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py
    
    """

    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.
    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.
    Returns
    -------
        A matplotlib ellipse artist
    """
    x = ps["PC-1"]
    y = ps["PC-2"]
    # print('PS-1')
    # print(x)
    # print('PS-2')
    # print(y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])

    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    theta = np.linspace(0, 2 * np.pi, size)

    ellipse_coords = np.column_stack([ell_radius_x * np.cos(theta), ell_radius_y * np.sin(theta)])

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    x_scale = np.sqrt(cov[0, 0]) * n_std
    x_mean = np.mean(x)

    # calculating the stdandard deviation of y ...
    y_scale = np.sqrt(cov[1, 1]) * n_std
    y_mean = np.mean(y)

    translation_matrix = np.tile([x_mean, y_mean], (ellipse_coords.shape[0], 1))
    rotation_matrix = np.array([[np.cos(np.pi / 4), np.sin(np.pi / 4)],
                                [-np.sin(np.pi / 4), np.cos(np.pi / 4)]])
    scale_matrix = np.array([[x_scale, 0],
                             [0, y_scale]])
    ellipse_coords = ellipse_coords.dot(rotation_matrix).dot(scale_matrix) + translation_matrix
    # print(ellipse_coords,"!!!")

    path = f'M {ellipse_coords[0, 0]}, {ellipse_coords[0, 1]}'
    for k in range(1, len(ellipse_coords)):
        path += f'L{ellipse_coords[k, 0]}, {ellipse_coords[k, 1]}'
    path += ' Z'
    return path
