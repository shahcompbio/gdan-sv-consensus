from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import argparse

def parse_args():
    desc = "Parse SV-supporting qnames (read names) from sniffles vcf"
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('-i', '--vcf', help="input sniffles vcf", required=True)
    p.add_argument('-s', '--sample', help="common sample ID", required=True)
    p.add_argument('-o', '--outdir', help="output dir for report", required=True)
    p.add_argument('--sources', help="list of SV sources", required=True, nargs='+')
    args = p.parse_args()
    return args

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

def qc_source_count(sources):
    n_sources = len(sources)
    if n_sources <= 1:
        raise ValueError(f'number of sources should be >1, not {n_sources}')
    elif n_sources > 3:
        raise ValueError(f'number of sources should not be >3, not {n_sources}')
        
def get_svs_from_union_vcf(vcf_path, sources):
    """From SURVIVOR union VCF, get regions specific to A, B, and regions in both A, and B
    """
    qc_source_count(sources)
    svtypes = ['DEL', 'DUP', 'INV', 'TRA']
    membership_count = {s:Counter() for s in svtypes}
    for line in open(vcf_path, 'r'):
        gts = {}
        members = {s: [0] * len(sources) for s in svtypes}
        if line.startswith('##'): continue
        elif line.startswith('#'):
            header = line.rstrip().split('\t')
            continue
        field = line.rstrip().split('\t')
        row = dict(zip(header, field))
        formats = row['FORMAT'].split(':')
        
        svtype_set = set()
        for ix, source in enumerate(sources):
            gts[source] = row[source].split(':')
            gts[source] = dict(zip(formats, gts[source]))
            svtype = gts[source]['TY'].split(',')[0]
            region = gts[source]['CO']
            flag_sv_found = int(svtype != 'NaN')
            #if flag_sv_found: svtype_set.add(svtype)
            if flag_sv_found:
                members[svtype][ix] = flag_sv_found
        for svtype in svtypes:
            membership = ''.join([str(m) for m in members[svtype]])
            membership_count[svtype][membership] += 1

    return dict(membership_count)

def get_set_name_map(sources):
    n_sources = len(sources)
    if n_sources == 2:
        key_list = [f'{bin(i)[2:]:0>2s}' for i in range(1, 2**n_sources)]
    elif n_sources == 3:
        key_list = [f'{bin(i)[2:]:0>3s}' for i in range(1, 2**n_sources)]
    else: raise ValueError(f'sources length should be 2 or 3, not {sources}')
    key_map = {}
    for key in key_list:
        set_name = [sources[i] if key[i]=='1' else '' for i in range(n_sources)]
        set_name = [v for v in set_name if v != '']
        key_map[key] = ' & '.join(set_name)
    key_list = sorted(key_map.values(), key=lambda x: len(x))
    key_rev_map = {v:k for (k, v) in key_map.items()}
    return key_rev_map, key_list

def plot_venn(args, membership_count):
    fig, axes = plt.subplots(2, 2, figsize=(9, 9))
    fig.suptitle(args.sample)
    svtypes = ['DEL', 'DUP', 'INV', 'TRA']
    cmap = {'DEL':'tab:blue', 'DUP':'tab:red', 'INV':'tab:purple', 'TRA':'grey'}
    for ix, svtype in enumerate(svtypes):
        yix = ix // 2
        xix = ix % 2
        ax = axes[yix][xix]
        ax.set_title(svtype, color=cmap[svtype])
        svtype_set_count = membership_count[svtype]
        if len(args.sources) == 2:
            venn2(subsets = svtype_set_count, set_labels=args.sources, ax=ax)
        elif len(args.sources) == 3:
            venn3(subsets = svtype_set_count, set_labels=args.sources, ax=ax)
        else:
            raise ValueError(f'args.sources={args.sources}')
    out_tag = f'{args.sample}.SV_union' 
    png_path = f'{args.outdir}/{out_tag}.venn.png'
    fig.savefig(png_path)

def write_report(args, membership_count):
    key_rev_map, key_list = get_set_name_map(args.sources)
    header_field = ['sample', 'svtype'] + key_list
    svtypes = ['DEL', 'DUP', 'INV', 'TRA']
    header = '\t'.join(header_field)
    out_tag = f'{args.sample}.SV_union' 
    with open(f'{args.outdir}/{out_tag}.report', 'w') as report:
        report.write(header + '\n')
        for svtype in svtypes:
            field = [args.sample, svtype]
            for key in key_list:
                name = key_rev_map[key]
                field.append(str(membership_count[svtype][name]))
            line = '\t'.join(field)
            report.write(line + '\n')


if __name__ == "__main__":
    args = parse_args()
    membership_count = get_svs_from_union_vcf(args.vcf, args.sources)
    write_report(args, membership_count)
    plot_venn(args, membership_count)


