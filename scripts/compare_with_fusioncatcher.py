import os
import glob
import argparse
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib_venn import venn2 

def parse_args():
    desc = "Compare fusioncatcher results with consensus bedpe"
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('--bedpe_dir', help="input bedpe dir", required=True)
    p.add_argument('--fusioncatcher_dir', help="input fusioncatcher dir", required=True)
    p.add_argument('-g', '--genes', nargs='+', help="genes to compare", required=True)
    p.add_argument('-s', '--samplesfile', help="samples file path")
    p.add_argument('-o', '--outdir', help="png output dir")
    args = p.parse_args()
    return args

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
    def __init__(self, sample, bedpe_dir='', fusioncatcher_dir=''):
        self.sample = sample
        self.short_sample = '-'.join(sample.split('-')[:2]) # CTSP-ACY0-TTP1-A-1-1-R-A78Z-41 -> CTSP-ACY0
        self.bedpe_dir = bedpe_dir
        self.fusioncatcher_dir = fusioncatcher_dir
        self.fusioncatcher = get_fusioncatcher(self.sample, indir=self.fusioncatcher_dir)
        self.bedpe_consensus = get_bedpe_consensus(self.short_sample, indir=self.bedpe_dir)

    def both_data_exist(self):
        flag = type(self.fusioncatcher) == pd.DataFrame and type(self.bedpe_consensus) == pd.DataFrame
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

def plot_overlap(sets, labels, title=None):
    set1, set2 = sets
    assert type(set1) == set, set1
    assert type(set2) == set, set2
    label1, label2 = labels
    fig, ax = plt.subplots(figsize=(4,4))
    venn2(subsets=(set1, set2), set_labels=(label1, label2), ax=ax)
    if title:
        ax.set_title(title)
    return fig, ax

def main():
    args = parse_args()
    bedpe_dir = args.bedpe_dir
    fs_dir = args.fusioncatcher_dir
    samples = [s.strip() for s in open(args.samplesfile)]
    samples_with_tra = defaultdict(list)
    samples_with_fusion = defaultdict(list)
    genes = args.genes

    for sample in samples:
        fc = FusionCatcher(sample, bedpe_dir=bedpe_dir, fusioncatcher_dir=fs_dir)
        if fc.both_data_exist():
            # n_exist += 1
            for gene in genes:
                if fc.fusion_for_gene_exists(gene):
                    samples_with_fusion[gene].append(sample)
                if fc.tra_for_gene_exists(gene):
                    samples_with_tra[gene].append(sample)

    with matplotlib.rc_context({'font.family':'Arial'}):
        for gene in genes:
            set1 = set(samples_with_tra[gene])
            set2 = set(samples_with_fusion[gene])
            label1 = 'Consensus SVs'
            label2 = 'FusionCatcher'
            fig, ax = plot_overlap((set1, set2), (label1, label2), title=gene)
            fig.savefig(f'{args.outdir}/fusion_vs_tra.{gene}.png')

if __name__ == '__main__':
    main()