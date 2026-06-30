import os, re, sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import wilcoxon, mannwhitneyu
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import multipletests
import argparse
# import param from param.py
from .param import CommandParams


def decide_sig(fdr, log2fc, log2fc_cutoff, fdr_cutoff):
    if log2fc > log2fc_cutoff and fdr < fdr_cutoff:
        return 'S-U'
    elif log2fc < -1 * log2fc_cutoff and fdr < fdr_cutoff:
        return 'S-D'
    elif log2fc > 0 and fdr < fdr_cutoff:
        return 'U'
    elif log2fc < 0 and fdr < fdr_cutoff:
        return 'D'


def compare_two_groups(df, a, b, method, fdr_cutoff=0.01, log2fc_cutoff=1,
                       max_miss_ratio_global=0.5, max_miss_ratio_group=0, min_sample_size=4):
    """
    Parameters
    ----------
    df: a data frame of expression matrix, index is gene/protein/id, column is sample IDs
    a: IDs of group 1
    b: IDs of group 2
    fdr_cutoff: FDR threshold for significant differential expression
    log2fc_cutoff: log2 fold change threshold for significant differential expression
    method: hypothesis method.
        'T-test(Unpaired)': Calculate the T-test for the means of two independent samples of scores.
        'T-test(Paired)':Calculate the T-test on two paired samples of scores.
        'Wilcoxon(Unpaired)': Perform the Mann-Whitney U rank test on two independent samples.
        'Wilcoxon(Paired)': Perform the Wilcoxon signed-rank test. It is a non-parametric version of the paired T-test.
    max_miss_ratio_global: maximum ratio of missing values among all samples(default: 0.5)
    max_miss_ratio_group: maximum ratio of missing values for each group (default: 1)
    min_sample_size: minimum samples size for comparison (default: 4)
    Returns
    -------
    rdf : a data frame of hypothesis test results with adjusted p-values
    """
    # check missing values and remove features having too many missing values
    rows = []
    methods = {
        'T-test(Unpaired)': ttest_ind,
        'T-test(Paired)': ttest_rel,
        'Wilcoxon(Unpaired)': mannwhitneyu,
        'Wilcoxon(Paired)': wilcoxon
    }
    # print(df.shape)
    for index, row in tqdm(df.iterrows()):

        miss_ratio = len([i for i in list(row) if pd.isna(i)]) / df.shape[1]
        if miss_ratio > max_miss_ratio_global:
            continue
        x = np.array(row[a], dtype=float)
        y = np.array(row[b], dtype=float)
        # print(x,y)
        if method in ['T-test(Paired)', 'Wilcoxon(Paired)']:
            # x = [x[i] for i in range(len(x)) if pd.notna(x[i]) and pd.notna(y[i])]
            # y = [y[i] for i in range(len(x)) if pd.notna(x[i]) and pd.notna(y[i])]
            z = [(i,j) for i,j in zip(x,y) if pd.notna(i) and pd.notna(j)]
            x = [i[0] for i in z]
            y = [i[1] for i in z]
            miss_ratio = (len(a) - len(x)) / len(a)
            # if miss_ratio > max_miss_ratio_group:
            #     continue
            if len(x) < min_sample_size:
                continue
            log2fc_median = np.median([i-j for i,j in zip(x,y)])
            log2fc_mean = np.mean([i-j for i,j in zip(x,y)])
        else:
            x = np.array([i for i in x if pd.notna(i)], dtype=float)
            miss_ratio_x = (len(a) - len(x)) / len(a)
            y = np.array([i for i in y if pd.notna(i)], dtype=float)
            miss_ratio_y = (len(b) - len(y)) / len(b)

            if miss_ratio_x > max_miss_ratio_group or miss_ratio_y > max_miss_ratio_group:
                # print(index, miss_ratio_x,miss_ratio_y,max_miss_ratio_group)
                continue
            # print(len(x),len(y))
            if len(x) < min_sample_size or len(y) < min_sample_size:
                continue
            log2fc_median = np.nanmedian(x) - np.nanmedian(y)
            log2fc_mean = np.nanmean(x) - np.nanmean(y)
        # print(index,len(x),len(y))

        test_result = methods[method](x, y)

        new_row = [index, log2fc_median, log2fc_mean, test_result[0], test_result[1]]
        rows.append(new_row)

    new_df = pd.DataFrame(rows, columns=['Feature', 'Log2FC(median)', 'Log2FC(mean)',
                                         '{}(Stats)'.format(method), '{}(P-value)'.format(method)])
    multiple_tests_result = multipletests(new_df['{}(P-value)'.format(method)], method='fdr_bh')
    new_df['FDR'] = multiple_tests_result[1]
    new_df['-Log10(FDR)'] = [-np.log10(i) for i in new_df['FDR']]

    new_df['Significance'] = [decide_sig(
        row['FDR'], row['Log2FC(median)'],log2fc_cutoff, fdr_cutoff
    ) for index, row in new_df.iterrows()]
    new_df = new_df.set_index('Feature')
    return new_df



if __name__ == "__main__":
    # run diff.py standalone
    # create argparse
    parser = argparse.ArgumentParser(description='Compare two groups of samples')
    # add input matrix
    parser.add_argument('--matrix', '-m', type=str, required=True,
                        help='expression matrix path')
    parser.add_argument('--param', '-p', type=str, required=True,
                        help='parameter file path')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='output directory')
    
    args = parser.parse_args()

    # get param path
    param_path = args.param
    matrix_path = args.matrix

    if not os.path.exists(param_path):
        print('Parameter file does not exist for diff.py')
        sys.exit()
    
    # generate param object
    params = CommandParams()
    params.read_from_config_path(param_path)

    # get expression matrix
    # run_diff(matrix_path,params)

    

    
