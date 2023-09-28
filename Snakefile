import os
import sys
import glob
import subprocess

import pandas as pd
import pybedtools
import wgs_analysis.refgenome as refgenome

if not os.path.exists(config['log_dir']): subprocess.run(f'mkdir -p {config["log_dir"]}', shell=True)
if not os.path.exists(config['tmp_dir']): subprocess.run(f'mkdir -p {config["tmp_dir"]}', shell=True)

def _get_samples(config, sources, sample_str_len=16):
    samples_per_source = {}
    for source in sources:
        if source == 'Broad':
            meta = pd.read_table(config['metadata'][source])
            samples = set(meta.iloc[:, 0].str.slice(0, sample_str_len))
        elif source == 'MSK':
            meta = pd.read_table(config['metadata'][source])
            samples = set(meta.iloc[:, 1].str.slice(0, sample_str_len))
        samples_per_source[source] = samples
    samples_intersection = samples_per_source[sources[0]]
    for source in sources[1:]:
        samples_intersection &= samples_per_source[source]
    return list(samples_intersection)

SOURCES = ['Broad', 'MSK', 'NYGC']
SAMPLES = ['CTSP-AD18-TTP1-A'] #_get_samples(config, SOURCES)

wildcard_constraints:
    source = '|'.join(SOURCES)

rule all:
    input:
        expand("results/{sample}/{sample}.SV_union.report", sample=SAMPLES),
        expand("results/{sample}/{sample}.SV_union.venn.png", sample=SAMPLES),
        expand("results/{sample}/{sample}.SV_intersection.vcf", sample=SAMPLES),
        expand("results/{sample}/{sample}.SV_union.vcf", sample=SAMPLES),
        expand('results/{sample}/{sample}.vcf_list.txt', sample=SAMPLES),
        expand('results/{sample}/{sample}.{source}.vcf', sample=SAMPLES, source=SOURCES),
        expand('results/{sample}/{sample}.{source}.bedpe', sample=SAMPLES, source=SOURCES),
        
def get_breakpoint_genes(gtf, brk, fai_path):
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
        return np.nan
    bed = pybedtools.BedTool(f'{chrom} {pos-1} {pos}', from_string=True)
    if strand == '+':
        closest = bed.closest(gtf, g=fai_path, 
                              fu=True, # get upstream
                              D="a") # report distance from "a[bed]"
    elif strand == '-':
        closest = bed.closest(gtf, g=fai_path, 
                              fd=True, # get upstream
                              D="a") # report distance from "a[bed]"
    else:
        raise ValueError(f'strand={strand} for {brk}')
    gene_info = [g[-2] for g in closest]
    gene_name = _extract_gene_name(brk, gene_info)
    return gene_name

def _extract_gene_name(brk, gene_info):
    # gene_info[0]
    # 'gene_id "ENSG00000227160.3"; gene_type "transcribed_unitary_pseudogene"; gene_name "THEM7P"; level 2; hgnc_id "HGNC:50386"; havana_gene "OTTHUMG00000154120.3";'
    gene_names = [g.split('gene_name "')[1].split('";')[0] for g in gene_info]
    gene_names = list(set(gene_names))
    gix = 0
    if len(gene_names) > 1:
        print(f'{brk}: more than one gene_names:{gene_names}', file=sys.stderr)
        gene_names_clean = [g for g in gene_names if g.count('.') == 0 and g.count('-') == 0]
        if len(gene_names_clean) >= 1:
            gix = gene_names.index(gene_names_clean[0])
    elif len(gene_names) == 0:
        return np.nan
    gene_name = gene_names[gix]
    return gene_name

def add_gene_names_to_svs(svs, gtf, fai_path):
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
        for ix, brk in enumerate(brks):
            ix += 1
            gene_name = get_breakpoint_genes(gtf, brk, fai_path)
            gene_names[ix].append(gene_name)

    for ix in [1, 2]:
        svs[f'gene{ix}'] = gene_names[ix]
    return svs

def _get_msk_svs_path(wildcards):
    paths = pd.read_table(config['metadata']['MSK'])
    paths = paths[paths['isabl_sample_id'].str.startswith(wildcards.sample)]
    paths = paths[paths['result_type']=='consensus_calls']
    if paths.shape[0] == 1:
        path = paths['result_filepath'].iloc[0]
        return path
    if paths.shape[0] == 0:
        return None
    else:
        raise ValueError(f'paths={paths} for sample={wildcards.sample}')

def get_msk_svs(path):
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
               'name', 'score', 'strand1', 'strand2', 'type',]
    if not path: # None
        return pd.DataFrame(columns=sv_cols)
    svs = pd.read_csv(path)
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
    svs['name'] = 'MSK' + svs['name']
    for sv_col in sv_cols:
        if sv_col not in svs.columns:
            svs[sv_col] = '.'
    chrom_same = svs['#chrom1'] == svs['chrom2']
    pos_reversed = svs['start1'] > svs['start2']
    rix = chrom_same & pos_reversed
    brk1_cols = ['#chrom1', 'start1', 'end1', 'strand1']
    brk2_cols = ['chrom2', 'start2', 'end2', 'strand2']
    svs.loc[rix, brk1_cols + brk2_cols] = svs.loc[rix, brk2_cols + brk1_cols].values
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

def _get_nygc_svs_path(wildcards):
    sample2path = {}
    nygc_sv_dir = config['metadata']['NYGC']
    paths = glob.glob(f'{nygc_sv_dir}/{wildcards.sample}*.sv.annotated.v7.somatic.high_confidence.final.bedpe')
    if len(paths) == 1:
        path = paths[0]
        return path
    elif len(paths) == 0:
        return None
    else:
        raise ValueError(f'paths={paths} for sample={wildcards.sample}')

def get_nygc_svs(svs_path):
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 
               'name', 'score', 'strand1', 'strand2', 'type', ]
    svs = pd.read_table(svs_path)
    svs.rename(columns={'#chr1':'#chrom1', 'chr2':'chrom2'}, inplace=True)
    svs['name'] = svs['tumor--normal'].str.slice(0, 16) + ':' + (svs.index + 1).astype(str)
    for sv_col in sv_cols:
        if sv_col not in svs.columns:
            svs[sv_col] = '.'
    return svs[sv_cols]

def _get_sv_type(row):
    svtype = row['type']
    chrom1, chrom2 = row['#chrom1'], row['chrom2']
    pos1, pos2 = row['start1'], row['end1']
    coord2 = f'{chrom2}:{pos2}'
    strands = f"{row.strand1}{row.strand2}"
    if svtype == 'TRA':
        if strands == '++':
            alt = f'N]{coord2}]'
        elif strands == '--':
            alt = f'[{coord2}[N'
        elif strands == '+-':
            alt = f'N[{coord2}['
        elif strands == '-+':
            alt = f']{coord2}]N'
        else:
            raise ValueError(f'strands={strands}')
        return alt
    else:
        return f'<{svtype}>'

def _get_vcf_meta_and_header(source, genome_version='hg38'):
    # genome_version = 'hg38'
    refgenome.set_genome_version(genome_version)
    chrom_infos = []
    for chrom, length in refgenome.info.chromosome_lengths.items():
        if genome_version == 'hg38':
            if 'chr' not in chrom: chrom = 'chr'+chrom
        chrom_info = '##contig=<ID={},length={}>'.format(chrom, length)
        chrom_infos.append(chrom_info)
    chrom_meta = '\n'.join(chrom_infos)

    meta_text = f'''##fileformat=VCFv4.2
{chrom_meta}
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend; Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">'''
    header_field = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', source]
    header = '\t'.join(header_field)
    return meta_text, header

def write_vcf_from_bedpe(sv_path, out_vcf, source='MSK'):
    with open(out_vcf, 'w') as out:
        meta_text, header = _get_vcf_meta_and_header(source=source, 
                                                     genome_version=config['ref']['genome_version'])
        out.write(meta_text + '\n')
        out.write(header + '\n')
        svs = pd.read_table(sv_path)
        svs[['start1', 'start2']] = svs[['start1', 'start2']].astype(int)

        for rix, row in svs.iterrows():
            chrom1, chrom2 = row['#chrom1'], row['chrom2']
            start1, start2 = row['start1'], row['start2']
            strand1, strand2 = row['strand1'], row['strand2']
            gene1, gene2 = row['gene1'], row['gene2']
            svname, svtype = row['name'], row['type']
            ref = 'N'
            alt = _get_sv_type(row)
            svid = f'{svname}__{gene1}__{gene2}'
            if (chrom1 == chrom2):
                if start2 < start1:
                    start1, start2 = start2, start1
                    strand1, strand2 = strand2, strand1
            chrom = chrom1
            pos = start1
            endpos = start2
            svlen = abs(endpos - pos)
            if row['type'] == 'DEL':
                svlen = -svlen
            strands = strand1 + strand2
            info = f'IMPRECISE;SVTYPE={svtype};END={endpos};STRAND={strands};CHR2={chrom2}'
            info += f';GENE1={gene1};GENE2={gene2};SOURCE={source}'
            qual = 60
            flt = 'PASS'
            fmt = 'GT:DV'
            gt = './.:99'
            field = [chrom, pos, svid, ref, alt, qual, flt, info, fmt, gt]
            field = [str(_) for _ in field]
            out.write('\t'.join(field) + '\n')

rule make_msk_bedpe:
    input:
        svs = _get_msk_svs_path,
    output:
        bedpe = 'results/{sample}/{sample}.MSK.bedpe',
    run:
        svs = get_msk_svs(input.svs)
        gtf = pybedtools.BedTool(config['ref']['gtf'])
        svs = add_gene_names_to_svs(svs, gtf, config['ref']['fai'])
        svs.to_csv(output.bedpe, sep='\t', index=False)
        
rule make_nygc_bedpe:
    input:
        svs = _get_nygc_svs_path,
    output:
        bedpe = 'results/{sample}/{sample}.NYGC.bedpe',
    run:
        svs = get_nygc_svs(input.svs)
        gtf = pybedtools.BedTool(config['ref']['gtf'])
        svs = add_gene_names_to_svs(svs, gtf, config['ref']['fai'])
        svs.to_csv(output.bedpe, sep='\t', index=False)
        
rule make_broad_bedpe:
    input:
        svs = config['metadata']['Broad'],
    output:
        bedpe = 'results/{sample}/{sample}.Broad.bedpe',
    run:
        svs = get_broad_svs(input.svs, wildcards.sample)
        gtf = pybedtools.BedTool(config['ref']['gtf'])
        svs = add_gene_names_to_svs(svs, gtf, config['ref']['fai'])
        svs.to_csv(output.bedpe, sep='\t', index=False)

rule make_vcf_from_bedpe:
    input:
        bedpe = 'results/{sample}/{sample}.{source}.bedpe',
    output:
        vcf = 'results/{sample}/{sample}.{source}.vcf',
    run:
        write_vcf_from_bedpe(input.bedpe, output.vcf, source=wildcards.source)
        
rule make_vcf_list:
    input:
        vcfs = expand('results/{{sample}}/{{sample}}.{source}.vcf', source=SOURCES),
    output:
        vcflist = 'results/{sample}/{sample}.vcf_list.txt',
    run:
        with open(output.vcflist, 'w') as out:
            for vcf_path in input.vcfs:
                out.write(vcf_path + '\n')
        
rule run_survivor:
    input:
        vcflist = 'results/{sample}/{sample}.vcf_list.txt',
    output:
        ivcf = "results/{sample}/{sample}.SV_intersection.vcf",
        uvcf = "results/{sample}/{sample}.SV_union.vcf",
    singularity: config['image']['survivor'],
    params:
        n_and = len(SOURCES),
        n_or = 1,
        diff = 30,
    shell:
        """
        /SURVIVOR/Debug/SURVIVOR merge {input} 200 {params.n_and} 1 1 0 {params.diff} {output.ivcf} &&
        /SURVIVOR/Debug/SURVIVOR merge {input} 200 {params.n_or}  1 1 0 {params.diff} {output.uvcf}
        """

rule summarize_survivor:
    input:
        uvcf = "results/{sample}/{sample}.SV_union.vcf",
    output:
        report = "results/{sample}/{sample}.SV_union.report",
        venn = "results/{sample}/{sample}.SV_union.venn.png",
    params:
        outdir = "results/{sample}",
        sources = ' '.join(SOURCES),
    shell:
        'python scripts/parse_survivor.py '
        '-i {input.uvcf} -o {params.outdir} '
        '-s {wildcards.sample} --sources {params.sources}'
