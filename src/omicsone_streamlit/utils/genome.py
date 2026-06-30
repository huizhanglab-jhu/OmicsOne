import pandas as pd
from tqdm import tqdm
import re
import numpy as np

def get_cytoband(chrom, pos, cytobands):
    """根据染色体和基因位置匹配染色体带"""
    band = cytobands[(cytobands["chrom"] == f"chr{chrom}") & 
                     (cytobands["start"] <= pos) & 
                     (cytobands["end"] >= pos)]
    return band["band"].iloc[0] if not band.empty else None

def add_band_info(gene_map, cytobands):
    for gene in tqdm(gene_map):
        chr = gene_map[gene]["chr"]
        pos = int(gene_map[gene]["offset"])
        band = get_cytoband(chr,pos,cytobands)
        gene_map[gene]['band'] = band
    return gene_map

def get_cytoband_map(cytoband_path):
    cytobands = pd.read_csv(cytoband_path,sep="\t", header=None)
    cytobands.columns = ['chrom','start','end','band','stain']
    cytoband_d = dict()
    for index,row in cytobands.iterrows():
        chrom = row['chrom']
        band = row['band']
        if pd.isna(band):
            continue
        key = f'{chrom[3:]}{band[0]}'
        if re.search('chr[\dXY]+$',chrom):
            start = row['start']
            end = row['end']
            if key not in cytoband_d:
                cytoband_d[key] = [start,end]
            else:
                if start < cytoband_d[key][0]:
                    cytoband_d[key][0] = start
                if end > cytoband_d[key][1]:
                    cytoband_d[key][1] = end
    return cytoband_d, cytobands

def calc_gistic(gistic_data,genes, samples, gene_map, start_map):
    
    genes = [i for i in genes if i.split(".")[0] in gene_map]
    
    gistic = gistic_data.loc[genes, samples]
    gistic.index = [i.split(".")[0] for i in gistic.index]
    
    gistic['chr'] = [gene_map.get(i)["chr"] for i in gistic.index]
    gistic['arm'] = [gene_map.get(i)["band"] for i in gistic.index]
    
    # print(gistic.head(2))
    gistic = gistic.drop_duplicates()
    gistic['chr.arm'] = gistic.apply(lambda row: row['chr'] + "." + row['arm'], axis=1)
    
    amplifications = (gistic.iloc[:, :-3] >= 1).sum(axis=1).reset_index(name='Count')
    deletions = (gistic.iloc[:, :-3] <= -1).sum(axis=1).reset_index(name='Count')
    
    # print(amplifications.head(2))
    # print(deletions.head(2))
    
    rows = []
    for index,row in amplifications.iterrows():
        gene = str(row['index'])
        chr = gene_map.get(gene)["chr"]
        start = start_map.get(chr)
        pos = start + gene_map.get(gene)["offset"]
        row['pos'] = pos
        row['chr'] = chr
        rows.append(row)

    amplifications_df = pd.DataFrame(rows)
    # amplifications_df

    rows = []
    for index,row in deletions.iterrows():
        gene = str(row['index'])
        chr = gene_map.get(gene)["chr"]
        start = start_map.get(chr)
        pos = start + gene_map.get(gene)["offset"]
        row['pos'] = pos
        row['chr'] = chr
        rows.append(row)

    deletions_df = pd.DataFrame(rows)
    
    lines_amp = []
    lines_del = []

    for index,row in amplifications_df.iterrows():
        x1 = row['pos']
        x2 = row['pos']
        y1 = 0
        y2 = row['Count']
        line = {"x1": x1, "y1": y1, "x2": x2, "y2": y2}
        lines_amp.append(line)

    for index,row in deletions_df.iterrows():
        x1 = row['pos']
        x2 = row['pos']
        y1 = 0
        y2 = row['Count'] * (-1)
        line = {"x1": x1, "y1": y1, "x2": x2, "y2": y2}
        lines_del.append(line)
        
    # Convert lines to a list of (x, y) pairs
    segments_amp = [([line["x1"], line["x2"]], [line["y1"], line["y2"]]) for line in lines_amp]

    # Prepare segments for LineCollection
    line_segments_amp = [((line["x1"], line["y1"]), (line["x2"], line["y2"])) for line in lines_amp]



    # Convert lines to a list of (x, y) pairs
    segments_del = [([line["x1"], line["x2"]], [line["y1"], line["y2"]]) for line in lines_del]

    # Prepare segments for LineCollection
    line_segments_del = [((line["x1"], line["y1"]), (line["x2"], line["y2"])) for line in lines_del]
    
    return line_segments_amp, line_segments_del
    
    

def cis_trans_corr_count(corr_df:pd.DataFrame, start_map: dict, cytoband_d:dict, corr_with : str)-> pd.DataFrame: 
    rows_pos_cis = []
    rows_pos_trans = []
    rows_neg_cis = []
    rows_neg_trans = []
    
    # print(corr_df.head(2))

    for index,row in tqdm(corr_df.iterrows()):
        # print(row['cnv.cytoband'], type(row['cnv.cytoband']))
        # print(row['rna.cytoband'], type(row['rna.cytoband']))
        # print(row['cnv.chr'], type(row['cnv.chr']))
        cnv_band = row['cnv.cytoband'][0]
        
        rna_band = None
        if corr_with == "RNA":
            rna_band = row['rna.cytoband'][0]
        elif corr_with == "Protein":
            rna_band = row['protein.cytoband'][0]
            
            
        cnv_chr = row['cnv.chr']
        rna_chr = None
        if corr_with == "RNA":
            rna_chr = row['rna.chr'] 
        elif corr_with == "Protein":
            rna_chr = row['protein.chr']
        
        cnv_gene = row['cnv.gene']
        
        rna_gene = None
        if corr_with == "RNA":
            rna_gene = row['rna.gene']
        elif corr_with == "Protein":
            rna_gene = row['protein.gene']
     
        corr = row['Correlation']
        if cnv_gene == rna_gene:
            continue
        if cnv_band == rna_band and cnv_chr == rna_chr:
            if corr > 0:
                rows_pos_cis.append(row)
            else:
                rows_neg_cis.append(row)
        else:
            if corr > 0:
                rows_pos_trans.append(row)
            else:
                rows_neg_trans.append(row)

    corr_pos_cis = pd.DataFrame(rows_pos_cis)
    corr_pos_trans = pd.DataFrame(rows_pos_trans)
    corr_neg_cis = pd.DataFrame(rows_neg_cis)
    corr_neg_trans = pd.DataFrame(rows_neg_trans)
    
    bands = []
    for i in list(range(1,23)) + ['X','Y']:
        for j in ['p','q']:
            bands.append(f'{i}{j}')
    
    if corr_pos_cis.shape[0] > 0:
        # corr_pos_cis_count = corr_pos_cis.groupby('cnv.band').size().reset_index(name='Count')
        corr_pos_cis_count = dict(corr_pos_cis.groupby('cnv.band').size())
    else:
        corr_pos_cis_count = dict()

    if corr_pos_trans.shape[0] > 0:
        corr_pos_trans_count = dict(corr_pos_trans.groupby('cnv.band').size())
    else:
        corr_pos_trans_count = dict()

    if corr_neg_cis.shape[0] > 0:
        corr_neg_cis_count = dict(corr_neg_cis.groupby('cnv.band').size())
    else:
        corr_neg_cis_count = dict()

    if corr_neg_trans.shape[0] > 0:
        corr_neg_trans_count = dict(corr_neg_trans.groupby('cnv.band').size())
    else:
        corr_neg_trans_count = dict()
        
    rows = []
    for arm in bands:
        # get count
        pos_cis_count = corr_pos_cis_count.get(arm,0)
        pos_trans_count = corr_pos_trans_count.get(arm,0)
        neg_cis_count = corr_neg_cis_count.get(arm,0)
        neg_trans_count = corr_neg_trans_count.get(arm,0)
        # get offset
        offset = start_map.get(arm[:-1])
        start,end = cytoband_d.get(arm)
        rows.append([arm, pos_cis_count, pos_trans_count, neg_cis_count, neg_trans_count, start + offset, end + offset])

    final_corr_count = pd.DataFrame(rows,
            columns = ['cnv.arm','pos.cis.count','pos.trans.count','neg.cis.count','neg.trans.count','cnv.arm.start','cnv.arm.end'])
    
    final_corr_count['count.sum'] = final_corr_count.apply(lambda row: np.sum([int(row['pos.cis.count']),int(row['pos.trans.count']),
                                                                         int(row['neg.cis.count']),int(row['neg.trans.count'])]), axis=1)
    return final_corr_count
        
    
    