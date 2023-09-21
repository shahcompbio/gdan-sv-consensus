import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse

def parse_args():
    desc = "Parse SV-supporting qnames (read names) from sniffles vcf"
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('-i', '--vcf', help="input sniffles vcf", required=True)
    p.add_argument('-s', '--sample', help="common sample ID", required=True)
    p.add_argument('-sa', '--sample_a', help="sample ID A", required=True)
    p.add_argument('-sb', '--sample_b', help="sample ID B", required=True)
    p.add_argument('-o', '--outdir', help="output dir for report", required=True)
    args = p.parse_args()
    return args

def get_svs_from_union_vcf(vcf_path, sample_a='Broad', sample_b='MSK'):
    """From SURVIVOR union VCF, get regions specific to A, B, and regions in both A, and B
    """
    onlyas = [] # ONT
    onlybs = [] # WGS
    bothas = [] # ONT
    bothbs = [] # WGS
    for line in open(vcf_path, 'r'):
        if line.startswith('##'): continue
        elif line.startswith('#'):
            header = line.rstrip().split('\t')
            continue
        field = line.rstrip().split('\t')
        row = dict(zip(header, field))
        formats = row['FORMAT'].split(':')
        
        gt_a = row[sample_a].split(':') # ONT
        gt_a = dict(zip(formats, gt_a))
        type_a = gt_a['TY']
        region_a = gt_a['CO']
        sv_a = f'{region_a}:{type_a}'
        
        gt_b = row[sample_b].split(':') # WGS
        gt_b = dict(zip(formats, gt_b))
        type_b = gt_b['TY']
        region_b = gt_b['CO']
        sv_b = f'{region_b}:{type_b}'
        
        if type_a != 'NaN' and type_b != 'NaN':
            bothas.append(sv_a)
            bothbs.append(sv_b)
        elif type_a == 'NaN':
            onlybs.append(sv_b)
        elif type_b == 'NaN':
            onlyas.append(sv_a)
        else:
            print(f'ERROR: id_a={id_a}, id_b={id_b}')
    return np.array(onlyas), np.array(onlybs), np.array(bothas), np.array(bothbs)

def process_sv_string(svstring, margin=1000):
    """Process SV string to IGV input form
    """
    chrompos, svtype = svstring.split(':')
    chrompos1, chrompos2 = chrompos.split('-')
    chrom1, pos1 = chrompos1.split('_')
    chrom2, pos2 = chrompos2.split('_')
    pos1 = int(pos1)
    pos2 = int(pos2)
    pos1_s, pos1_e = pos1-margin, pos1+margin
    pos2_s, pos2_e = pos2-margin, pos2+margin
    return f'{chrom1}:{pos1_s}-{pos1_e} {chrom2}:{pos2_s}-{pos2_e} {svtype}'


if __name__ == "__main__":
    args = parse_args()

    onlyas, onlybs, bothas, bothbs = get_svs_from_union_vcf(args.vcf, args.sample_a, args.sample_b)
    assert len(bothas) == len(bothbs)
    out_tag = f'{sample}.SV_union' # CTSP-AD1H-TTP1-A.SV_union.venn.png

    with open(f'{args.outdir}/{out_tag}.report', 'w') as report:
        nA_not_B = len(onlyas) 
        nB_not_A = len(onlybs)
        nA_and_B = len(bothas)
        nA = nA_not_B + nA_and_B
        nB = nB_not_A + nA_and_B
        nA_or_B = nA_not_B + nB_not_A + nA_and_B
        header = ['A', 'B', 'A-B', 'B-A', 'A&B', 'A|B']
        field = [nA, nB, nA_not_B, nB_not_A, nA_and_B, nA_or_B]
        field = [str(x) for x in field]
        report.write('\t'.join(header) + '\n')
        report.write('\t'.join(field) + '\n')

        # assume A: ONT, B: WGS
        plt.figure(figsize=(4,4))
        plt.title(args.sample)
        png_path = f'{args.outdir}/{out_tag}.venn.png'
        venn2(subsets = (nA_not_B, nB_not_A, nA_and_B), set_labels=('ONT', 'WGS'))
        plt.savefig(png_path)

