import os
import sys
import glob
import subprocess

import tabix
import numpy as np
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

CHROMS = ['chr'+str(c) for c in range(1, 22+1)] + ['chrX', 'chrY']
SOURCES = ['Broad', 'MSK', 'NYGC']
SAMPLES = ['CTSP-AD18-TTP1-A'] #_get_samples(config, SOURCES)

wildcard_constraints:
    source = '|'.join(SOURCES)

rule all:
    input:
        #expand('results/{sample}/{sample}.consensus_stringent.bedpe', sample=SAMPLES),
        #expand('results/{sample}/{sample}.consensus_lenient.bedpe', sample=STRINGENT),
        'results/gtf/protein_coding.gtf.gz',
        'results/gtf/protein_coding.gtf.gz.tbi',
        expand("results/{sample}/{sample}.SV_union.report", sample=SAMPLES),
        expand("results/{sample}/{sample}.SV_union.venn.png", sample=SAMPLES),
        expand("results/{sample}/{sample}.SV_union.vcf", sample=SAMPLES),
        expand('results/{sample}/{sample}.vcf_list.txt', sample=SAMPLES),
        expand('results/{sample}/{sample}.{source}.vcf', sample=SAMPLES, source=SOURCES),
        expand('results/{sample}/{sample}.{source}.bedpe', sample=SAMPLES, source=SOURCES),
        
rule grep_and_sort_gtf:
    input:
        gtf = config['ref']['gtf'],
    output:
        chrom_gtf = temp(os.path.join(config['tmp_dir'], '{chrom}.gtf')),
    shell: 
        'zcat {input.gtf} | grep -P "^{wildcards.chrom}\t" | '
        'sort -k4,5 -n > {output.chrom_gtf}'

rule merge_gtf:
    input:
        expand(os.path.join(config['tmp_dir'], '{chrom}.gtf'), chrom=CHROMS),
    output:
        gtf = 'results/gtf/genes.gtf'
    shell:
        'cat {input} > {output.gtf}'

rule filter_gtf:
    input: 
        gtf = 'results/gtf/genes.gtf'
    output: 
        gtf = 'results/gtf/protein_coding.gtf.gz',
        tbi = 'results/gtf/protein_coding.gtf.gz.tbi',
    shell:
        'cat {input.gtf} | grep protein_coding | bgzip -c > {output.gtf} && '
        'tabix -f -p gff {output.gtf}'

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

def _get_sv_type(row):
    svtype = row['type']
    chrom1, chrom2 = row['#chrom1'], row['chrom2']
    pos1, pos2 = row['end1'], row['end2']
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

def _get_sv_path(wildcards):
    if wildcards.source == 'MSK':
        path = _get_msk_svs_path(wildcards)
    elif wildcards.source == 'NYGC':
        path = _get_nygc_svs_path(wildcards)
    elif wildcards.source == 'Broad':
        path = config['metadata']['Broad']
    else:
        raise ValueError(f'source={source}')
    return path

rule make_bedpe:
    input:
        svs = _get_sv_path,
        gtf = 'results/gtf/protein_coding.gtf.gz',
    output:
        bedpe = 'results/{sample}/{sample}.{source}.bedpe',
    params:
        fai = config['ref']['fai'],
    shell:
        'python scripts/make_bedpe_from_svs.py -i {input.svs} -o {output.bedpe} '
        '-g {input.gtf} -f {params.fai} -s {wildcards.source} --sample {wildcards.sample}'
        
#rule make_msk_bedpe:
#    input:
#        svs = _get_msk_svs_path,
#        gtf = 'results/gtf/protein_coding.gtf.gz',
#    output:
#        bedpe = 'results/{sample}/{sample}.MSK.bedpe',
#    run:
#        svs = get_msk_svs(input.svs)
#        gtf = pybedtools.BedTool(input.gtf)
#        tbx = tabix.open(input.gtf)
#        svs = add_gene_names_to_svs(svs, gtf, tbx, config['ref']['fai'])
#        svs.to_csv(output.bedpe, sep='\t', index=False)
#        
#rule make_nygc_bedpe:
#    input:
#        svs = _get_nygc_svs_path,
#        gtf = 'results/gtf/protein_coding.gtf.gz',
#    output:
#        bedpe = 'results/{sample}/{sample}.NYGC.bedpe',
#    run:
#        svs = get_nygc_svs(input.svs)
#        gtf = pybedtools.BedTool(input.gtf)
#        tbx = tabix.open(input.gtf)
#        svs = add_gene_names_to_svs(svs, gtf, tbx, config['ref']['fai'])
#        svs.to_csv(output.bedpe, sep='\t', index=False)
#        
#rule make_broad_bedpe:
#    input:
#        svs = config['metadata']['Broad'],
#        gtf = 'results/gtf/protein_coding.gtf.gz',
#    output:
#        bedpe = 'results/{sample}/{sample}.Broad.bedpe',
#    run:
#        svs = get_broad_svs(input.svs, wildcards.sample)
#        gtf = pybedtools.BedTool(input.gtf)
#        tbx = tabix.open(input.gtf)
#        svs = add_gene_names_to_svs(svs, gtf, tbx, config['ref']['fai'])
#        svs.to_csv(output.bedpe, sep='\t', index=False)

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
        vcf = "results/{sample}/{sample}.SV_union.vcf",
    singularity: config['image']['survivor'],
    params:
        n_and = len(SOURCES),
        n_or = 1,
        diff = 30,
    shell:
        """
        /SURVIVOR/Debug/SURVIVOR merge {input} 200 {params.n_or}  1 1 0 {params.diff} {output.vcf}
        """

rule summarize_survivor:
    input:
        vcf = "results/{sample}/{sample}.SV_union.vcf",
    output:
        report = "results/{sample}/{sample}.SV_union.report",
        venn = "results/{sample}/{sample}.SV_union.venn.png",
    params:
        outdir = "results/{sample}",
        sources = ' '.join(SOURCES),
    shell:
        'python scripts/parse_survivor.py '
        '-i {input.vcf} -o {params.outdir} '
        '-s {wildcards.sample} --sources {params.sources}'

rule create_consensus_bedpe:
    input:
        vcf = "results/{sample}/{sample}.SV_union.vcf",
    output:
        stringent = 'results/{sample}/{sample}.consensus_stringent.bedpe',
        lenient = 'results/{sample}/{sample}.consensus_lenient.bedpe',
    shell:
        'python scripts/create_consensus_bedpe.py '
        '-i {input.vcf} -os {output.stringent} -ol {output.lenient}'
