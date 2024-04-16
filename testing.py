import os
import sys
import yaml
from unittest import TestCase

import tabix
import pybedtools
import pandas as pd

os.chdir('/rtsess01/compute/juno/shah/users/chois7/projects/GDAN/DLBCL/pipeline')
sys.path.append('scripts')
config = yaml.load(open('config.yaml'), Loader=yaml.Loader)
from make_bedpe_from_svs import get_overlapping_genes, get_downstream_gene, get_upstream_gene, add_gene_names_to_svs
from make_consensus_bedpe import make_consensus_bedpe


gtf_path = 'results/gtf/protein_coding.gtf.gz'
fai_path = config['ref']['fai']
assert(os.path.exists(gtf_path)), gtf_path

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