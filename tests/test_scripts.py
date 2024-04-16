import os
import sys
import yaml
from unittest import TestCase

import tabix
import pybedtools
import pandas as pd

sys.path.append('scripts')
config = yaml.load(open('config.yaml'), Loader=yaml.Loader)
from make_bedpe_from_svs import get_overlapping_genes, get_downstream_gene, get_upstream_gene, add_gene_names_to_svs
from make_consensus_bedpe import make_consensus_bedpe

gtf_path = 'results/gtf/protein_coding.gtf.gz'
fai_path = config['ref']['fai']
assert(os.path.exists(gtf_path)), gtf_path

class TestGeneFetching(TestCase):
    def test_get_overlapping_genes(self):
        tbx = tabix.open(gtf_path)
        brk = ('chr14', 105863215, '+')
        self.assertTrue(get_overlapping_genes(tbx, brk) == 'IGHJ6')

    def test_get_downstream_gene(self):
        gtf = pybedtools.BedTool(gtf_path)
        brks = [('chr14', 105862938, '-'), ('chr14', 105872938, '-')]
        self.assertTrue(get_downstream_gene(gtf, brks, fai_path) == 'IGHJ6')

    def test_get_upstream_gene(self):
        gtf = pybedtools.BedTool(gtf_path)
        brks = [('chr14', 105862938, '-'), ('chr14', 105872938, '-')]
        self.assertTrue(get_upstream_gene(gtf, brks, fai_path) == 'IGHD7-27')
    
    def test_add_gene_names_to_svs(self):
        gtf = pybedtools.BedTool(gtf_path)
        tbx = tabix.open(gtf_path)
        svs = pd.DataFrame({'#chrom1': {0: 'chr8'},
                            'start1': {0: 127736662},
                            'end1': {0: 127736663},
                            'chrom2': {0: 'chr14'},
                            'start2': {0: 105862937},
                            'end2': {0: 105862938},
                            'name': {0: 'TEST'},
                            'score': {0: 1},
                            'strand1': {0: '+'},
                            'strand2': {0: '-'},
                            'type': {0: 'TRA'}})
        svs = add_gene_names_to_svs(svs, gtf, tbx, fai_path)
        self.assertTrue(svs['gene1'].values[0] == 'MYC')
        self.assertTrue(svs['gene2'].values[0] == 'IGHJ6')

class TestConsensus(TestCase):
    def test_make_consensus_bedpe(self):
        vcf_path = 'tests/test_data.vcf'
        sources = ['Broad', 'MSK', 'NYGC']
        goi = {'BCL2', 'BCL6', 'MYC'}
        vcf_cols = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + sources
        df = pd.read_table(vcf_path, names=vcf_cols, comment='#')
        bedpe = make_consensus_bedpe(df, sources, goi, consensus_cutoff=2)
        self.assertTrue(bedpe.iloc[:2]['score'].tolist() == [3, 2])
        self.assertTrue(bedpe.shape[0] == 3)
        self.assertTrue(bedpe.iloc[-1]['gene1'] == 'MYC')
        self.assertTrue(bedpe.iloc[-1]['gene2'] == 'IGHJ6')
        self.assertTrue(bedpe.iloc[-1]['score'] == 1) # salvage successful