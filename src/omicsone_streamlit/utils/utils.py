import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np

import re, os, sys
# from .mylogger import mylogger

import matplotlib.ticker as ticker
from tqdm import tqdm

# imports
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from sklearn.model_selection import train_test_split
from statsmodels.stats.multitest import multipletests
from scipy import stats

# -*- coding: utf-8 -*-

'''
Created on 2011-8-24
@author: JerryKwan, Yingwei Hu

'''

import logging
import logging.handlers

import os

LEVELS = {'NOSET': logging.NOTSET,
          'DEBUG': logging.DEBUG,
          'INFO': logging.INFO,
          'WARNING': logging.WARNING,
          'ERROR': logging.ERROR,
          'CRITICAL': logging.CRITICAL}


class LogMgr:
    def __init__(self, logpath, markpath):
        self.LOG = logging.getLogger('log')
        loghdlr1 = logging.handlers.RotatingFileHandler(logpath, "a", 0, 1)
        fmt1 = logging.Formatter("%(asctime)s %(threadName)-10s %(message)s", "%Y-%m-%d %H:%M:%S")
        loghdlr1.setFormatter(fmt1)
        self.LOG.addHandler(loghdlr1)
        self.LOG.setLevel(logging.INFO)

        self.MARK = logging.getLogger('mark')
        loghdlr2 = logging.handlers.RotatingFileHandler(markpath, "a", 0, 1)
        fmt2 = logging.Formatter("%(message)s")
        loghdlr2.setFormatter(fmt2)
        self.MARK.addHandler(loghdlr2)
        self.MARK.setLevel(logging.INFO)

    def error(self, msg):
        if self.LOG is not None:
            self.LOG.error(msg)

    def info(self, msg):
        if self.LOG is not None:
            self.LOG.info(msg)

    def debug(self, msg):
        if self.LOG is not None:
            self.LOG.debug(msg)

    def mark(self, msg):
        if self.MARK is not None:
            self.MARK.info(msg)

#
# mylogger = LogMgr("mylog", "mymark")


def print_result(func):
    def wrapper(*args, **kw):
        result = func(*args, **kw)
        # args_str = str(args)
        kw_str = str(kw)
        status = 'FAIL' if result is None else 'SUCCESS'
        # message = 'call {0}():{1},{2}'.format(func.__name__, kw_str, status)
        message = "call {0}".format(func.__name__)
        print(message)
        return result

    return wrapper


def save_table(df, out, **kwargs):
    # index = kwargs.get('index', True)
    bn, ext = os.path.splitext(out)
    if ext == ".csv":
        df.to_csv(out, sep=",", **kwargs)
    elif ext == ".tsv":
        df.to_csv(out, sep="\t", **kwargs)
    elif ext == ".xlsx":
        df.to_excel(out, **kwargs)


def read_txt(filename):
    if not os.path.exists(filename):
        # mylogger.error('ERR:Failed to read file: {0}, File Not Found ').format(filename)
        return None
    try:
        with open(filename, 'r') as f:
            x = f.readlines()
            x = [i for i in x if not re.match("#", i)]
            x = [re.sub("\n$", '', i) for i in x]
            x = [str(i).strip() for i in x]
            x = [i for i in x if i != ""]
            f.close()
            if len(x) == 0:
                # mylogger.error('WARN:Failed to read file: {0}, File does not contain valid information ').format(
                #     filename)
                return None
            else:
                return x
    except Exception as e:
        # mylogger.error('ERR:Failed to read file: {0}, because {1} ').format(filename, str(e))
        return None


def read_xlsx(filename, sheet_name=0):
    if not os.path.exists(filename):
        return None
    try:
        df = pd.read_excel(filename, sheet_name=sheet_name)
        return df
    except Exception as e:
        # mylogger.error('ERR:Failed to read file: {0}, because {1} '.format(filename, str(e)))
        return None


def read_tsv(filename, sep="\t"):
    if not os.path.exists(filename):
        return None
    try:
        df = pd.read_csv(filename, sep=sep)
        return df
    except Exception as e:
        # mylogger.error('ERR:Failed to read file: {0}, because {1} ').format(filename, str(e))
        return None


def read_csv(filename, sep=","):
    return read_tsv(filename, sep=sep)


def read_table(filename):
    # print("filename",filename)
    bn, ext = os.path.splitext(filename)
    # print(bn,ext)
    read_table_methods = {
        ".txt": read_txt,
        ".csv": read_csv,
        ".tsv": read_tsv,
        ".xlsx": read_xlsx
    }
    return read_table_methods[ext](filename)


# def myboxplot(df, gene, dpi=300):
#     tumor = df[df.index == gene].loc[:, cols['malignant']].T
#     tumor['State'] = 'Tumor'
#     print('tumor', tumor[pd.notnull(tumor[gene])].shape[0])
#     normal = df[df.index == gene].loc[:, cols['nonmalignant']].T
#     normal['State'] = 'Non-tumor'
#     print('normal', normal[pd.notnull(normal[gene])].shape[0])
#     x = pd.concat([tumor, normal])
#     sns.set_style("white")
#     fig = plt.figure(figsize=(2, 3), dpi=dpi)
#     sns.boxplot(y=gene, x='State', data=x, palette="Greys_r")
#     ax = sns.swarmplot(y=gene, x='State', data=x, color='1', edgecolor="black", linewidth=0.5)
#     ax.set_ylabel('log2 ratio value', fontsize=16)
#     ax.set_xlabel('', fontsize=16)
#     plt.tick_params(labelsize=12)
#     t, p = stats.ttest_ind(list([i for i in tumor[gene] if not np.isnan(i)]),
#                            list([i for i in normal[gene] if not np.isnan(i)]))
#     if p < 0.01:
#         ax.set_title('{0}\nP-value:{1:.2e}'.format(gene, p), fontsize=18)
#     else:
#         ax.set_title('{0}\nP-value:{1:.2f}'.format(gene, p), fontsize=18)


def myboxplot2(df, gene, groups):
    g1 = df[df.index == gene].loc[:, groups[0]].T
    g1['State'] = 'ST1'
    g2 = df[df.index == gene].loc[:, groups[1]].T
    g2['State'] = 'ST2'
    g3 = df[df.index == gene].loc[:, groups[2]].T
    g3['State'] = 'ST3'
    x = pd.concat([g1, g2, g3])
    sns.set_style("white")
    fig = plt.figure(figsize=(6, 4), dpi=200)
    sns.boxplot(y=gene, x='State', data=x, palette="Greys_r")
    ax = sns.swarmplot(y=gene, x='State', data=x, color='1', edgecolor="black", linewidth=0.5)
    ax.set_ylabel('log2 ratio value', fontsize=20)
    ax.set_xlabel('Pathological Status', fontsize=20)
    plt.tick_params(labelsize=16)
    t12, p12 = stats.ttest_ind([i for i in list(g1[gene]) if pd.notnull(i)],
                               [i for i in list(g2[gene]) if pd.notnull(i)])

    t23, p23 = stats.ttest_ind([i for i in list(g2[gene]) if pd.notnull(i)],
                               [i for i in list(g3[gene]) if pd.notnull(i)])
    print(p12, p23)
    #     if p12 < 0.05:
    ax.set_title('{0}, P-value:{1:.2e},{2:.2e}'.format(gene, p12, p23), fontsize=20)


#     else:
#         ax.set_title('{0}, P-value:{1:.2f}'.format(gene,p),fontsize=20)


def box_plot(df, gene, groups, group_names, **kwargs):
    figure_title = kwargs.get('figure_title', 'box_plot')
    figure_title_font_size = int(kwargs.get('figure_title_font_size', 16))
    figure_width = int(kwargs.get('figure_width', 5))
    figure_height = int(kwargs.get('figure_height', 5))
    figure_dpi = int(kwargs.get('figure_dpi', 300))
    figure_x_label = kwargs.get('figure_x_label', 'Group')
    figure_x_label_font_size = int(kwargs.get('figure_x_label_font_size', 12))
    figure_y_label = kwargs.get('figure_y_label', 'Gene')
    figure_y_label_font_size = int(kwargs.get('figure_y_label_font_size', 12))
    figure_tick_font_size = int(kwargs.get('figure_tick_font_size', 10))
    fout = kwargs.get('fout', None)

    gs = []
    for i in range(len(groups)):
        g = df[df.index == gene].loc[:, groups[i]].T
        g['State'] = group_names[i]
        gs.append(g)
    x = pd.concat(gs)
    sns.set_style('white')
    fig = plt.figure(figsize=(figure_width, figure_height), dpi=figure_dpi)
    sns.boxplot(y=gene, x='State', data=x, palette="Greys_r")
    ax = sns.swarmplot(y=gene, x='State', data=x, color='1', edgecolor="black", linewidth=0.5)
    ax.set_ylabel(figure_y_label, fontsize=figure_y_label_font_size)
    ax.set_xlabel(figure_x_label, fontsize=figure_x_label_font_size)
    plt.tick_params(labelsize=figure_tick_font_size)
    plt.xticks(rotation=45)
    ax.set_title(figure_title, fontsize=figure_title_font_size)
    if fout is not None:
        plt.savefig(fout, bbox_inches='tight', dpi=figure_dpi)
    else:
        plt.show()


def mk_groups(orig_df, group_column, ordered_groups=None):
    df = orig_df.copy(deep=True)
    df = df.applymap(lambda i: 'NA' if pd.isnull(i) else i)
    if ordered_groups is None:
        group_names = df[group_column]
        group_names = sorted(list(set(group_names)))
        groups = [list(df[df[group_column] == i].index) for i in group_names]
        return groups, group_names
    else:
        groups = [list(df[df[group_column] == i].index) for i in ordered_groups]
        return groups, ordered_groups


def get_pos(x, name):
    y = [i for i in range(len(x)) if name == x[i]]
    if len(y) > 0:
        return y[0]
    else:
        return None


def volcano_plot(result, group1, group2,
                 x_col='log2FC', x_cutoff=0,
                 y_col='adjusted_p', y_cutoff=0.05, y_neglg10=True,
                 **fig_kwargs):
    sns.set_style("white")
    default_colors = dict(pos='r', neg='b')
    colors = fig_kwargs.get('colors', default_colors)
    fout = fig_kwargs.get('fout', None)

    y_col_new = '-log10({})'.format(y_col)
    if y_neglg10 and (y_col_new not in result.columns.values):
        result[y_col_new] = [-np.log10(i) for i in result[y_col]]

    result_pos = result[result[x_col] > x_cutoff]
    result_neg = result[result[x_col] < x_cutoff]

    figsize = fig_kwargs.get('figsize', (4, 4))
    dpi = fig_kwargs.get('dpi', 300)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    sns.regplot(x=x_col,
                y=y_col_new if y_neglg10 else y_col,
                data=result_pos, fit_reg=False, color=colors['pos'],
                scatter_kws={'s': 1, 'alpha': 0.3})
    sns.regplot(x=x_col,
                y=y_col_new if y_neglg10 else y_col,
                data=result_neg, fit_reg=False, color=colors['neg'],
                scatter_kws={'s': 1, 'alpha': 0.3})

    highlight = fig_kwargs.get('highlight', False)
    if highlight:
        sig_pos_x = fig_kwargs.get('sig_pos_fc', 2)
        sig_neg_x = fig_kwargs.get('sig_neg_fc', 2)
        sig_pos_y = fig_kwargs.get('sig_pos_p', 0.01)
        sig_neg_y = fig_kwargs.get('sig_neg_p', 0.01)
        result_sig_pos = result[(result[x_col] >= np.log2(sig_pos_x)) & (result[y_col] <= sig_pos_y)]
        result_sig_neg = result[(result[x_col] <= -np.log2(sig_neg_x)) & (result[y_col] <= sig_neg_y)]

        sns.regplot(x=x_col,
                    y=y_col_new if y_neglg10 else y_col,
                    data=result_sig_pos, fit_reg=False, color='r',
                    scatter_kws={'s': 1, 'alpha': 1})
        sns.regplot(x=x_col,
                    y=y_col_new if y_neglg10 else y_col,
                    data=result_sig_neg, fit_reg=False, color='b',
                    scatter_kws={'s': 1, 'alpha': 1})

        plt.axhline(y=-np.log10(sig_pos_y), color='lightgray', lw=1, ls='-.', label='P-value<0.05')
        plt.axvline(x=-np.log2(sig_neg_x), color='lightgray', lw=1, ls='-.', label='P-value<0.05')
        plt.axvline(x=np.log2(sig_pos_x), color='lightgray', lw=1, ls='-.', label='P-value<0.05')

    g1 = fig_kwargs.get('group_name1', 'Group1')
    g2 = fig_kwargs.get('group_name2', 'Group2')
    method = fig_kwargs.get('method', "")
    plt.ylabel('{} test(-log10(FDR))'.format(method))
    plt.xlabel('Median difference(log2FC, {0} vs {1})'.format(g1, g2))

    plt.title('{0} vs {1}: {3} vs {4}'.format(g1, g2, result.shape[0], len(group1), len(group2)))
    x_max = np.max([abs(i) for i in result[x_col]]) + 1
    x_min = fig_kwargs.get('x_min', -1 * x_max)
    x_max = fig_kwargs.get('x_max', x_max)
    plt.xlim(x_min, x_max)
    # plt.show()
    if fout is not None:
        plt.savefig(fout, bbox_inches='tight', dpi=300)
    else:
        plt.show()


def bar_plot(df, **kwargs):
    figure_width = int(kwargs.get('figure_width', 5))
    figure_height = int(kwargs.get('figure_height', 5))
    figure_dpi = int(kwargs.get('figure_dpi', 300))
    figure_color = kwargs.get('figure_color', 'red')
    fig = plt.figure(figsize=(figure_width, figure_height), dpi=figure_dpi)
    sns.barplot(x=df.index, y='Count', data=df, color=figure_color)
    plt.xticks(rotation='vertical')
    fout = kwargs.get('fout', None)
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        if ext not in ['.png', 'jpg']:
            fout = bn + '.png'
        plt.savefig(fout, bbox_inches='tight', dpi=figure_dpi)
    else:
        plt.show()


def heat_map(df, mask=None, **kwargs):
    # Set up the matplotlib figure
    figure_width = int(kwargs.get('figure_width', 11))
    figure_height = int(kwargs.get('figure_height', 9))
    f, ax = plt.subplots(figsize=(figure_width, figure_height))

    # Generate a custom diverging colormap
    # cmap = sns.diverging_palette(220, 10, as_cmap=True)
    cmap = sns.diverging_palette(240, 10, n=9, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(df, mask=mask, cmap=cmap, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    fout = kwargs.get('fout', None)
    figure_dpi = int(kwargs.get('figure_dpi', 300))
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        if not re.search('heatmap$', bn):
            bn += "_heatmap"
        if ext not in ['png', 'jpg']:
            fout = bn + '.png'
        plt.savefig(fout, bbox_inches='tight', dpi=figure_dpi)
    else:
        plt.show()


def dist_plot(df, **kwargs):
    # Set up the matplotlib figure
    figure_width = int(kwargs.get('figure_width', 11))
    figure_height = int(kwargs.get('figure_height', 9))
    f, ax = plt.subplots(figsize=(figure_width, figure_height))

    sns.distplot(df['Value'])
    fout = kwargs.get('fout', None)
    figure_dpi = int(kwargs.get('figure_dpi', 300))
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        if not re.search('distplot$', bn):
            bn += "_distplot"
        if ext not in ['png', 'jpg']:
            fout = bn + '.png'
        plt.savefig(fout, bbox_inches='tight', dpi=figure_dpi)
    else:
        plt.show()


def pair_plot(df, **kwargs):
    # Set up the matplotlib figure
    figure_width = int(kwargs.get('figure_width', 11))
    figure_height = int(kwargs.get('figure_height', 9))
    f, ax = plt.subplots(figsize=(figure_width, figure_height))

    sns.pairplot(df, kind='reg', plot_kws=dict(scatter_kws={"s": 0.5}))
    fout = kwargs.get('fout', None)
    figure_dpi = int(kwargs.get('figure_dpi', 300))
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        if not re.search('pairplot$', bn):
            bn += "_pairplot"
        if ext not in ['png', 'jpg']:
            fout = bn + '.png'
        plt.savefig(fout, bbox_inches='tight', dpi=figure_dpi)
    else:
        plt.show()


def vectorize_stage(i, remove_substage=True):
    if remove_substage:
        i = re.sub('[ABC]', '', i)
    d = {'I': 1, 'II': 2, 'III': 3, 'IV': 4}
    return d.get(i, '[NA]')


def vectorize_grade(i):
    i = re.sub('G', '', i)
    if re.match('[0-9]+$', i):
        return int(i)
    else:
        return '[NA]'


def vectorize_col_stage(col):
    d = {'I': 1, 'II': 2, 'III': 3, 'IV': 4}
    r = col.apply(lambda i: d[re.sub('[ABC]', '', i)])
    return r


def vectorize_col_grade(col):
    r = col.apply(lambda i: vectorize_grade(i))
    return r


def vectorize_col_default(col):
    values = sorted(list(set(col)))
    #     print(values)
    d = dict(zip(values, range(1, len(values) + 1)))
    #     print(d)
    r = col.apply(lambda i: d[i])
    return r


def vectorize_col_numerical(col):
    r = [np.nan if pd.isnull(i) else float(i) for i in col]
    return r


vectorization = {
    'stage': vectorize_col_stage,
    'grade': vectorize_col_grade,
    'numerical': vectorize_col_numerical
}

vectorization_types = {
    'tumor_stage': 'stage',
    'tumor_grade': 'grade',
    'age_at_diagnosis': 'numerical',
}


def vectorize(df):
    df2 = df.copy(deep=True)
    for col in df2.columns.values:
        if col in vectorization_types:
            method = vectorization_types[col]
            # print(method)
            df2[col] = vectorization[method](df2[col])
            # print(df2[col])
        else:
            df2[col] = vectorize_col_default(df2[col])
    return df2


def count_plot(df, x, **kwargs):
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    # sns.distplot(info_tumors_validCellularity['% Tumor Cellularity'],kde=False)
    plt.title(kwargs.get('figure_title', 'Count {}'.format(x)))

    sns.countplot(data=df.sort_values(x), x=x)
    temp = plt.ylabel(kwargs.get('ylable', 'Count'))
    plt.xlabel(kwargs.get('xlabel', ''))
    fout = kwargs.get('fout', None)
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        figure_path = bn + ".png"
        plt.savefig(figure_path, bbox_inches='tight', dpi=300)
    else:
        plt.show()


def read_groups(df, param):
    groups = param.get('groups', 'all')
    if groups == 'all':
        groups = list(df.columns.values)
    else:
        groups = re.split('/', groups)
    return groups


def save_table(df, param, tag=''):
    fout = param.get('fout', None)
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        table_path = bn + tag + ".csv"
        df.to_csv(table_path, index=True)


def save_figure(plt, param, **kwargs):
    fout = param.get('fout', None)
    tag = kwargs.get('tag', '')
    if fout is not None:
        bn, ext = os.path.splitext(fout)
        figure_path = bn + tag + ".png"
        plt.savefig(figure_path, bbox_inches='tight', **kwargs)
    else:
        plt.show()


def histogram(values, param):
    figure_size = param.get('figure_size', (5, 5))
    figure_dpi = param.get('figure_dpi', 300)
    fig, ax = plt.subplots(figsize=figure_size, dpi=figure_dpi)
    figure_title = param.get('figure_title', 'Histogram')

    plt.hist(values, bins=20, color='white', edgecolor='black')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    # print(plt.xticks())
    median = np.median(values)
    if median is not None:
        plt.axvline(float(median), color='b', linestyle='dashed', linewidth=2)
    figure_title += '\n{:.2f}'.format(median)
    plt.title(figure_title, fontsize=18)
    plt.xlabel('Coefficient of Variation (CV)', fontsize=16)
    plt.ylabel('Number', fontsize=16)
    plt.xlim(0, 1)
    kwargs = {'dpi': figure_dpi}
    save_figure(plt, param, **kwargs)
    # if out is not None:
    #     plt.savefig(out, bbox_inches='tight', **kwargs)


def read_gmt(path):
    lines = []
    with open(path, 'r') as f:
        lines = f.readlines()
        f.close()
    data = []
    for line in lines:
        values = re.split('\t', line)
        pathway_name = values[0]
        pathway_genes = values[2:]
        data.append([pathway_name, len(pathway_genes), ','.join(pathway_genes)])
    pathways_df = pd.DataFrame(data, columns=['Name', 'GeneNumber', 'Genes'])
    return pathways_df


def check_pathway_significance(orig_df, orig_info, param, enrichr_res, tag=""):
    shared_samples = set(orig_df.columns.values) & set(orig_info.index)
    df = orig_df[shared_samples]
    df.to_csv("C:/Users/yhu39/temp/df.csv")
    info = orig_info[orig_info.index.isin(shared_samples)]
    group = param.get('phenotype_for_pathway', 'daystodeath_or_LFU')

    if group == 'daystodeath_or_LFU':
        info = info[info['vital_status'] == 'DECEASED']

    er = enrichr_res.sort_values('P-value')
    term = er['Term'].iloc[0]

    default_enrich_dir = "/home/yingwei/newdrive/bitbucket/omicsone/omicsone/res/databases/enrichment"
    enrich_dir = param.get('enrich_database', default_enrich_dir)
    gmt_path = enrich_dir + "/combined_mkhr.gmt"
    gmt = read_gmt(gmt_path)
    pathway = gmt[gmt['Name'] == term]
    genes = re.split(',', str(pathway['Genes'].iloc[0]).strip())
    genes = [i for i in genes if i in df.index]

    pg_df = df[df.index.isin(genes)]
    pg_df = pg_df[[i for i in pg_df.columns.values if i in info.index]]

    reg_res = []
    df2 = df.T
    # df2.to_csv("C:/Users/yhu39/temp/df2.csv")
    df2.columns = [re.sub('[^A-Z0-9a-z]', '_', i) for i in df2.columns.values]
    df2 = df2.applymap(float)
    info2 = info.merge(df2, left_index=True, right_index=True)
    # info2.to_csv("C:/Users/yhu39/temp/info2.csv")

    confounder = param.get('confounder', 'age_at_diagnosis')

    for i in tqdm(df2.columns.values):
        if confounder == "":
            formula = '{0} ~ {1}'.format(group, i)
            lm1 = smf.ols(formula=formula, data=info2).fit()
            i_coef = dict(lm1.params)[i]
            i_p = dict(lm1.pvalues)[i]
            reg_res.append([np.nan, i_coef, np.nan, i_p, i])

        else:
            formula = '{0} ~ {1} + {2}'.format(group, confounder, i)
            lm1 = smf.ols(formula=formula, data=info2).fit()
            confounder_coef = dict(lm1.params)[confounder]
            confounder_p = dict(lm1.pvalues)[confounder]
            i_coef = dict(lm1.params)[i]
            i_p = dict(lm1.pvalues)[i]
            reg_res.append([confounder_coef, i_coef, confounder_p, i_p, i])
    reg_df = pd.DataFrame(reg_res, columns=['{}_coef'.format(confounder), 'Gene_coef',
                                            '{}_p'.format(confounder), 'Gene_p', 'Index'])
    multiple_tests_method = 'fdr_bh'
    reg_df['Gene_p_adjust'] = multipletests(reg_df['Gene_p'], method=multiple_tests_method)[1]
    reg_df["neglogP"] = [(-1) * i for i in reg_df['Gene_p_adjust'].apply(np.log10)]

    reg_df = reg_df.set_index('Index')
    reg_df[term] = [term if i in genes else "Other Genes" for i in reg_df.index]

    save_table(reg_df, param, tag='_sig_pathway')
    a1 = reg_df[reg_df.index.isin(genes)]['neglogP']
    a2 = reg_df[~reg_df.index.isin(genes)]['neglogP']

    stest = stats.mannwhitneyu(a1, a2, alternative='greater')
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    sns.boxplot(x=term, y='neglogP', data=reg_df)
    plt.xlabel('Gene')
    plt.ylabel('-log10(P) of one-sided Mannwhitneyu Test ')
    plt.title('Association with {}, P-value={:.2e}'.format(group, list(stest)[1]))
    tag = "_{}_sig_pathway".format(tag) if tag != "" else "_sig_pathway"
    save_figure(plt, param, tag=tag)


def to_r_str(mylist):
    r_str = ",".join(["'{}'".format(i) for i in mylist])
    return r_str


#
# def quick_import_param(param_path):
#     from .param import OmicsOneParameter
#     from .preprocess import preprocess
#
#     if param_path is None or not os.path.exists(param_path):
#         print('Error! Fail to load parameters from {}'.format(param_path))
#         sys.exit(1)
#     param = OmicsOneParameter(param_path)
#
#     # check output directory
#     if not os.path.exists(param['report']['out_dir']):
#         wd = os.path.dirname(os.path.abspath(__file__))
#         out_dir = wd + "/../res/" + param['report']['out_dir']
#         if os.path.exists(out_dir):
#             param['report']['out_dir'] = out_dir
#         else:
#             # print(out_dir)
#             print('Output Directory Not Found: {0} or {1}'.format(
#                 param['report']['out_dir'], out_dir
#             ))
#             try:
#                 os.mkdir(out_dir)
#                 param['report']['out_dir'] = out_dir
#                 print('Create output directory: {}'.format(
#                     out_dir
#                 ))
#             except:
#                 print('Fail to create output directory: {}'.format(
#                     out_dir
#                 ))
#                 return None
#
#     # check file existence
#     for i in param['data']:
#         fn = param['data'][i]['path']
#         method = "FILE_EXISTENCE_CHECK"
#         if os.path.exists(fn):
#             result = "SUCCESS"
#             param['data'][i][method] = result
#             print('{0}: input={1},result={2}'.format(method, fn, result))
#         else:
#             wd = os.path.dirname(os.path.abspath(__file__))
#             fn = wd + "/../res/" + fn
#             param['data'][i]['path'] = fn
#             if os.path.exists(fn):
#                 result = "SUCCESS"
#                 param['data'][i][method] = result
#                 print('{0}: input={1},result={2}'.format(method, fn, result))
#             else:
#                 result = "FAIL"
#                 param['data'][i][method] = result
#                 print('{0}: input={1},result={2}'.format(method, fn, result))
#
#     # preprocess data
#     for data_id in param['data']:
#         if param['data'][data_id]['FILE_EXISTENCE_CHECK'] != "SUCCESS":
#             continue
#         # the file is not involved in the analysis
#         if data_id not in param['method']:
#             continue
#         if data_id in param['method']:
#             preprocess_param = param['method'][data_id].get('preprocess', None)
#             if preprocess_param:
#                 fn = param['data'][data_id]['path']
#                 df = preprocess(file_name=fn, param=preprocess_param)
#                 param['data'][data_id]['preprocessed'] = df
#
#     return param


def get_genes_in_pathway(pathway_name, **param):
    wd = os.path.dirname(os.path.abspath(__file__))
    default_gmt = wd + "/../res/databases/enrichment/combined_mkhr.gmt"
    gmt_path = param.get('gmt', default_gmt)
    gmt = read_gmt(gmt_path)
    pathway = gmt[gmt['Name'] == pathway_name]
    genes = re.split(',', str(pathway['Genes'].iloc[0]).strip())
    return genes


def path_full_split(path):
    folder, bn = os.path.split(path)
    bn, ext = os.path.splitext(bn)
    return folder, bn, ext

def standard_read_expression_file(path):
    folder, bn, ext = path_full_split(path)
    df = None
    if ext == ".xlsx":
        df = pd.read_excel(path, engine="openpyxl", index_col=0)
    elif ext == ".txt" or ext == ".tsv":
        df = pd.read_csv(path, sep="\t", index_col=0)
    elif ext == ".csv":
        df = pd.read_csv(path, index_col=0)
    return df
