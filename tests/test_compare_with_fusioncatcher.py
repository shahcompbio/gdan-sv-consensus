import sys
sys.path.append('scripts')

import glob

import pytest
import matplotlib.pyplot as plt
from matplotlib_venn import venn2 
import pandas as pd
from compare_with_fusioncatcher import (
    FusionCatcher, return_data_if_exists_else_none, get_fusioncatcher, get_bedpe_consensus,
    plot_overlap
)

def test_return_data_if_exists_else_none():
    paths = []
    assert return_data_if_exists_else_none([]) == None

def test_get_fusioncatcher():
    sample = 'CTSP-ACY0-TTP1-A-1-1-R-A78Z-41'
    indir = '../data/rnaseq/fusions/fusioncatcher'
    df = get_fusioncatcher(sample, indir)
    assert type(df) == pd.DataFrame
    assert df.shape[0] > 0, df

def test_get_bedpe_consensus():
    sample = 'CTSP-ACXZ-TTP1-A'
    indir='../data/wgs/breakpoints/bedpe_consensus'
    df = get_bedpe_consensus(sample, indir)
    assert type(df) == pd.DataFrame
    assert df.shape[0] > 0, df

class TestFusionCatcher:
    sample = 'CTSP-ACY0-TTP1-A-1-1-R-A78Z-41'
    fc = FusionCatcher(sample)

    def test_init(self):
        fc = self.fc
        assert hasattr(fc, 'fusioncatcher')
        assert hasattr(fc, 'bedpe_consensus')
        assert type(fc.both_data_exist()) == bool

    def test_gene_exists(self):
        fc = self.fc
        genes = ['MYC', 'BCL2', 'BCL2']
        for gene in genes:
            assert type(fc.tra_for_gene_exists(gene)) == bool
            assert type(fc.fusion_for_gene_exists(gene)) == bool

    def test_fusion_for_gene_exists(self):
        fc = self.fc
        assert fc.fusion_for_gene_exists('BCL2')

@pytest.mark.mpl_image_compare
def test_plot_overlap():
    sets, labels = [{'a', 'b'}, {'b', 'c'}], ['set1', 'set2']
    assert len(sets) == len(labels)
    fig, ax = plot_overlap(sets, labels)
    return fig