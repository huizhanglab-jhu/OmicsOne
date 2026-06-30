import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse



def plot_pca1(df, palette='Dark2',n=1,figsize=(18, 12),xlim=(-20,30),ylim=(-20,30),title='PCA Plot'):
    # np.random.seed(0)
    # df = pd.DataFrame({
    #     'PC1': np.random.randn(100),
    #     'PC2': np.random.randn(100),
    #     'PC3': np.random.randn(100),
    #     'PC4': np.random.randn(100),
    #     'Group': np.random.choice(range(1, 11), 100)
    # })
    
    # 假设最后一列是组别信息
    groups = list(df.iloc[:, -1])
    data = df.iloc[:, :-1]

    # 进行PCA分析
    pca = PCA(n_components=4)
    principalComponents = pca.fit_transform(data)

    df_pca = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2', 'PC3', 'PC4'])
    df_pca.index = data.index
    df_pca['Group'] = groups
    
    df_pca = df_pca.sort_values(by='Group')

    # 函数：绘制椭圆
    def plot_ellipse(x, y, ax,n_std=2.0, **kwargs):
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
    
    pc1_variance = pca.explained_variance_[0]  # 第一个特征值
    pc2_variance = pca.explained_variance_[1]  # 第一个特征值
    total_variance = sum(pca.explained_variance_)  # 总方差
    pc1_variance_ratio = pca.explained_variance_ratio_[0]  # PC1 解释的方差占比
    pc2_variance_ratio = pca.explained_variance_ratio_[1]  # PC1 解释的方差占比

    # 创建图形
    if n == 1:
        fig, ax = plt.subplots(figsize=figsize)
        sns.scatterplot(data=df_pca, x='PC1', y='PC2', hue='Group', palette=palette, 
                        ax=ax, s=50, alpha=0.7, legend='full')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # 绘制椭圆
        n_groups = len(df_pca['Group'].unique())
        for index, group in enumerate(sorted(list(df_pca['Group'].unique()))):
            group_data = df_pca[df_pca['Group'] == group]
            plot_ellipse(group_data['PC1'], group_data['PC2'], ax=ax, 
                        edgecolor=sns.color_palette(palette, n_groups)[index])

        ax.set_xlabel(f'PC-1({pc1_variance_ratio*100:.2f}%)')
        ax.set_ylabel(f'PC-2({pc2_variance_ratio*100:.2f}%)')
    else:
        fig, axes = plt.subplots(2, 3, figsize=figsize)

        pairs = [('PC1', 'PC2'), ('PC1', 'PC3'), ('PC1', 'PC4'), ('PC2', 'PC3'), ('PC2', 'PC4'), ('PC3', 'PC4')]
        for (ax, (pc_x, pc_y)) in zip(axes.flatten(), pairs):
            sns.scatterplot(data=df_pca, x=pc_x, y=pc_y, hue='Group', palette=palette, 
                            ax=ax, s=50, alpha=0.7, legend='full')

            # 绘制椭圆
            n_groups = len(df_pca['Group'].unique())
            for index, group in enumerate(sorted(list(df_pca['Group'].unique()))):
                group_data = df_pca[df_pca['Group'] == group]
                plot_ellipse(group_data[pc_x], group_data[pc_y], ax=ax, 
                            edgecolor=sns.color_palette(palette, n_groups)[index])

            ax.set_xlabel(pc_x)
            ax.set_ylabel(pc_y)
    plt.title(title)

    plt.tight_layout()
    # plt.show()
    
    return fig, pca



def plot_pca(df, palette='Dark2'):
    
    # 假设最后一列是组别信息
    groups = list(df.iloc[:, -1])
    data = df.iloc[:, :-1]

    # 进行PCA分析
    pca = PCA(n_components=4)
    principalComponents = pca.fit_transform(data)
    
    explained_var = pca.explained_variance_ratio_ * 100  # 转换为百分比

    df_pca = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2', 'PC3', 'PC4'])
    df_pca.index = data.index
    df_pca['Group'] = groups
    
    df_pca = df_pca.sort_values(by='Group')

    # 函数：绘制椭圆
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

    # 创建图形
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    pairs = [('PC1', 'PC2'), ('PC1', 'PC3'), ('PC1', 'PC4'), ('PC2', 'PC3'), ('PC2', 'PC4'), ('PC3', 'PC4')]
    for (ax, (pc_x, pc_y)) in zip(axes.flatten(), pairs):
        sns.scatterplot(data=df_pca, x=pc_x, y=pc_y, hue='Group', palette=palette, 
                        ax=ax, s=50, alpha=0.7, legend='full')

        # 绘制椭圆
        n_groups = len(df_pca['Group'].unique())
        for index, group in enumerate(sorted(list(df_pca['Group'].unique()))):
            group_data = df_pca[df_pca['Group'] == group]
            plot_ellipse(group_data[pc_x], group_data[pc_y], ax=ax, 
                         edgecolor=sns.color_palette(palette, n_groups)[index])

        
        # 提取 PC 的编号索引
        pc_x_index = int(pc_x.replace("PC", "")) - 1
        pc_y_index = int(pc_y.replace("PC", "")) - 1

        # 设置含方差比例的坐标标签
        ax.set_xlabel(f"{pc_x} ({explained_var[pc_x_index]:.1f}%)")
        ax.set_ylabel(f"{pc_y} ({explained_var[pc_y_index]:.1f}%)")

    plt.tight_layout()

    
    return fig, df_pca
