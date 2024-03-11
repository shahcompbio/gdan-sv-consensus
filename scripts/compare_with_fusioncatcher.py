import os
import glob
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2 

def return_data_if_exists_else_none(paths):
    exists = False
    data = pd.DataFrame()
    for path in paths:
        if os.path.exists(path):
            exists = True
            df = pd.read_table(path)
            data = pd.concat([data, df])
    if exists:
        return data
    return None

def get_fusioncatcher(sample, indir='../data/rnaseq/fusions/fusioncatcher'):
    short_sample = '-'.join(sample.split('-')[:2]) # CTSP-ACY0-TTP1-A-1-1-R-A78Z-41 -> CTSP-ACY0
    paths = glob.glob(f'{indir}/{short_sample}*_final-list_candidate-fusion-genes.txt')
    data = return_data_if_exists_else_none(paths)
    return data

def get_bedpe_consensus(sample, indir='../data/wgs/breakpoints/bedpe_consensus'):
    short_sample = '-'.join(sample.split('-')[:2]) # CTSP-ACY0-TTP1-A-1-1-R-A78Z-41 -> CTSP-ACY0
    paths = glob.glob(f'{indir}/{sample}*.SV_consensus.bedpe')
    data = return_data_if_exists_else_none(paths)
    return data

class FusionCatcher:

    def __init__(self, sample):
        self.sample = sample
        self.short_sample = '-'.join(sample.split('-')[:2]) # CTSP-ACY0-TTP1-A-1-1-R-A78Z-41 -> CTSP-ACY0
        self.fusioncatcher = get_fusioncatcher(self.sample)
        self.bedpe_consensus = get_bedpe_consensus(self.short_sample)

    def both_data_exist(self):
        flag = type(self.fusioncatcher.shape[0]) == pd.DataFrame and type(self.bedpe_consensus.shape[0]) == pd.DataFrame
        return flag

    def tra_for_gene_exists(self, gene):
        tra = self.bedpe_consensus[self.bedpe_consensus['type'] == 'TRA']
        flag = (gene in tra['gene1'].values) or (gene in tra['gene2'].values)
        return flag

    def fusion_for_gene_exists(self, gene):
        df = self.fusioncatcher
        flag = (gene in df['Gene_1_symbol(5end_fusion_partner)'].values) or (gene in df['Gene_2_symbol(3end_fusion_partner)'].values)
        return flag

def shorten_sample(sample):
    short_sample = '-'.join(sample.strip().split('-')[:2]) # CTSP-ACY0-TTP1-A-1-1-R-A78Z-41 -> CTSP-ACY0
    return short_sample

def plot_overlap(sets, labels):
    set1, set2 = sets
    label1, label2 = labels
    fig, ax = plt.subplots(figsize=(4,4))
    venn2(subsets=(set1, set2), set_labels=(label1, label2), ax=ax)
    return fig, ax