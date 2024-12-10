import sys
import argparse

import tabix
import pybedtools
import numpy as np
import pandas as pd


def parse_args():
    desc = "Parse SV-supporting qnames (read names) from sniffles vcf"
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('-i', '--sv', help="input sv file path", required=True)
    p.add_argument('-o', '--bedpe', help="output bedpe path", required=True)
    p.add_argument('-g', '--gtf', help="gtf path", required=True)
    p.add_argument('-f', '--fai', help="fai path", required=True)
    p.add_argument('-s', '--source', help="svs source", required=True)
    p.add_argument('--sample', help="sample ID")
    args = p.parse_args()
    return args

def get_nygc_svs(svs_path):
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
               'name', 'score', 'strand1', 'strand2', 'type', ]
    svs = pd.read_table(svs_path)
    svs.rename(columns={'#chr1':'#chrom1', 'chr2':'chrom2'}, inplace=True)
    svs['name'] = svs['tumor--normal'].str.slice(0, 16) + '_' + (svs.index + 1).astype(str)
    for sv_col in sv_cols:
        if sv_col not in svs.columns:
            svs[sv_col] = '.'
    return svs[sv_cols]

def get_msk_svs(path):
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
               'name', 'score', 'strand1', 'strand2', 'type',]
    if not path: # None
        return pd.DataFrame(columns=sv_cols)
    svs = pd.read_csv(path, dtype={'chromosome_1':str, 'chromosome_2':str, 
        'position_1':int, 'position_2':int, 'prediction_id':str})
    col_map = {
        'chromosome_1': '#chrom1', 'chromosome_2': 'chrom2',
        'position_1': 'end1', 'position_2': 'end2',
        'strand_1': 'strand1', 'strand_2': 'strand2', 'prediction_id': 'name',
    }
    type_map = {
        'duplication': 'DUP', 'deletion': 'DEL', 'inversion': 'INV',
        'translocation': 'TRA', 'insertion': 'INS'
    }
    svs.rename(columns=col_map, inplace=True)
    svs['start1'] = svs['end1'] - 1
    svs['start2'] = svs['end2'] - 1
    svs['type'] = svs['type'].map(type_map)
    svs['name'] = 'MSK' + svs['name'].astype(str)
    svs['score'] = svs['num_split'] + svs['num_unique_reads']
    for sv_col in sv_cols:
        if sv_col not in svs.columns:
            svs[sv_col] = '.'
    chrom_same = svs['#chrom1'] == svs['chrom2']
    pos_reversed = svs['start1'] > svs['start2']
    rix = chrom_same & pos_reversed
    brk1_cols = ['#chrom1', 'start1', 'end1', 'strand1']
    brk2_cols = ['chrom2', 'start2', 'end2', 'strand2']
    svs.loc[rix, brk1_cols + brk2_cols] = svs.loc[rix, brk2_cols + brk1_cols].values
    svs[['start1', 'start2', 'end1', 'end2']] = svs[['start1', 'start2', 'end1', 'end2']].astype(int)
    return svs[sv_cols]

def get_broad_svs(path, sample): # tumor_submitter_id      individual      chrom1  start1  end1    chrom2  start2  end2    sv_id   tumreads        strand1 strand2 svclass svmethod
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
               'name', 'score', 'strand1', 'strand2', 'type',]
    col_map = {'chrom1':'#chrom1', 'sv_id':'name', 'tumreads':'score', 'svclass':'type'}
    svs = pd.read_table(path, sep='\t')
    svs = svs[svs['tumor_submitter_id']==sample]
    svs.rename(columns=col_map, inplace=True)
    svs['type'] = svs['type'].replace({'h2hINV':'INV', 't2tINV':'INV'})
    svs['#chrom1'] = 'chr'+svs['#chrom1']
    svs['chrom2'] = 'chr'+svs['chrom2']
    return svs[sv_cols]

def _extract_gene_name(gene_info):
    # gene_info[0]
    # 'gene_id "ENSG00000227160.3"; gene_type "transcribed_unitary_pseudogene"; gene_name "THEM7P"; level 2; hgnc_id "HGNC:50386"; havana_gene "OTTHUMG00000154120.3";'
    gene_names = [g.split('gene_name "')[1].split('";')[0] for g in gene_info]
    gene_names = list(set(gene_names))
    gix = 0
    if len(gene_names) > 1:
        print(f'{gene_info}: more than one gene_names:{gene_names}', file=sys.stderr)
        gene_names_clean = [g for g in gene_names if g.count('.') == 0 and g.count('-') == 0]
        if len(gene_names_clean) >= 1:
            gix = gene_names.index(gene_names_clean[0])
    elif len(gene_names) == 0:
        return ''
    gene_name = gene_names[gix]
    return gene_name

def get_overlapping_genes(tbx, brk):
    """Get genes overlapping a breakpoint 

    :param tbx: PyTabix GTF
    :type tbx: pybedtools.bedtool.BedTool
    :param brk: (chrom, pos, strand) tuple
    :type brk: tuple
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    chrom, pos, _ = brk
    pos = int(pos)
    results = tbx.query(chrom, pos-1, pos)
    gene_info = [g[8] for g in results]
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def get_closest_genes(gtf, brk, fai_path):
    """Get genes near a breakpoint with strand taken into account

    :param gtf: BedTools GTF
    :type gtf: pybedtools.bedtool.BedTool
    :param brk: (chrom, pos, strand) tuple
    :type gtf: tuple
    :param fai_path: path to fai file
    :type fai_path: str
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    chrom, pos, strand = brk
    canon_chroms = set([f'chr{c}' for c in range(1,22+1)] + ['chrX', 'chrY'])
    if chrom not in canon_chroms: 
        return ''
    bed = pybedtools.BedTool(f'{chrom} {pos-1} {pos}', from_string=True)
    closest = bed.closest(gtf, g=fai_path,
                          fu=(strand=='+'), # get upstream if '+'
                          fd=(strand=='-'), # get downstream if '-'
                          D="a") # report distance from "a[bed]"=>brk        
    gene_info = [g[11] for g in closest]
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def get_upstream_gene(gtf, brks, fai_path):
    """Get genes upstream of a breakpoint that is inside an adjacency

    :param gtf: BedTools GTF
    :type gtf: pybedtools.bedtool.BedTool
    :param brks: 2-element list of (chrom, pos, strand) tuple
    :type gtf: list
    :param fai_path: path to fai file
    :type fai_path: str
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    brk1, brk2 = brks
    chrom1, pos1, _ = brk1 
    chrom2, pos2, _ = brk2 # look upstream to this breakpoint
    assert chrom1 == chrom2, (chrom1, chrom2)
    canon_chroms = set([f'chr{c}' for c in range(1,22+1)] + ['chrX', 'chrY'])
    if chrom2 not in canon_chroms: 
        return ''
    bed = pybedtools.BedTool(f'{chrom2} {pos2-1} {pos2}', from_string=True)
    closest = bed.closest(gtf, g=fai_path,
                          k=1, # look for one gene
                          fu=True, # get upstream
                          D="a") # report distance from "a[bed]"        
    gene_info = [g[11] for g in closest if pos1 <= pos2-abs(int(g[-1]))] # dist g[-1] is negative
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def get_downstream_gene(gtf, brks, fai_path):
    """Get genes downstream of a breakpoint that is inside an adjacency

    :param gtf: BedTools GTF
    :type gtf: pybedtools.bedtool.BedTool
    :param brks: 2-element list of (chrom, pos, strand) tuple
    :type gtf: list
    :param fai_path: path to fai file
    :type fai_path: str
    ...
    :return: gene_gtf tuple
    :rtype: tuple
    """
    brk1, brk2 = brks
    chrom1, pos1, _ = brk1 # look downstream to this breakpoint
    chrom2, pos2, _ = brk2 
    assert chrom1 == chrom2, (chrom1, chrom2)
    canon_chroms = set([f'chr{c}' for c in range(1,22+1)] + ['chrX', 'chrY'])
    if chrom1 not in canon_chroms: 
        return ''
    bed = pybedtools.BedTool(f'{chrom1} {pos1-1} {pos1}', from_string=True)
    closest = bed.closest(gtf, g=fai_path,
                          k=1, # look for one gene
                          fd=True, # get downstream
                          D="a") # report distance from "a[bed]"        
    gene_info = [g[11] for g in closest if pos1+abs(int(g[-1])) <= pos2] # g[-1] is positive
    gene_name = _extract_gene_name(gene_info)
    return gene_name

def add_gene_names_to_svs(svs, gtf, tbx, fai_path):
    gene_names = {1:[], 2:[]}
    svs = svs.copy()
    for rix, row in svs.iterrows():
        if rix % 10 == 0: print(f'progress: {rix}/{svs.shape[0]}')
        chrom1, chrom2 = row['#chrom1'], row['chrom2']
        pos1, pos2 = row['end1'], row['end2']
        strand1, strand2 = row['strand1'], row['strand2']
        svtype = row['type']
        assert (svtype=='TRA') or (svtype!='TRA' and chrom1==chrom2 and pos1<=pos2), row
        brks = [(chrom1, pos1, strand1), (chrom2, pos2, strand2)]
        if svtype in {'TRA', 'INV'}:
            for ix, brk in enumerate(brks):
                bix = ix + 1
                gene_name = get_overlapping_genes(tbx, brk)
                if len(gene_name) == 0: # ''
                    gene_name = get_closest_genes(gtf, brk, fai_path)
                gene_names[bix].append(gene_name)
        else:
            gene_name_1 = get_overlapping_genes(tbx, brks[0])
            if len(gene_name_1) == 0:
                gene_name_1 = get_downstream_gene(gtf, brks, fai_path)
            gene_names[1].append(gene_name_1)
            gene_name_2 = get_overlapping_genes(tbx, brks[1])
            if len(gene_name_2) == 0:
                gene_name_2 = get_upstream_gene(gtf, brks, fai_path)
            gene_names[2].append(gene_name_2)

    for ix in [1, 2]:
        svs[f'gene{ix}'] = gene_names[ix]
    return svs

def get_svs(path, source, sample=None):
    if source == 'MSK':
        svs = get_msk_svs(path)
    elif source == 'NYGC':
        svs = get_nygc_svs(path)
    elif source == 'Broad':
        svs = get_broad_svs(path, sample)
    else:
        raise ValueError(f'source={source}')
    return svs


if __name__ == "__main__":
    args = parse_args()
    svs = get_svs(args.sv, args.source, sample=args.sample)
    gtf = pybedtools.BedTool(args.gtf)
    tbx = tabix.open(args.gtf)
    svs = add_gene_names_to_svs(svs, gtf, tbx, args.fai)
    svs.to_csv(args.bedpe, sep='\t', index=False)

