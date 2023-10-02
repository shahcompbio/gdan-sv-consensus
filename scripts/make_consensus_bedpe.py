import re
import argparse
import pandas as pd

def parse_args():
    desc = "Parse SV-supporting qnames (read names) from sniffles vcf"
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('-i', '--vcf', help="input SURVIVOR vcf path", required=True)
    p.add_argument('-o', '--bedpe', help="output bedpe path", required=True)
    p.add_argument('-s', '--sources', nargs='+', help="svs sources", required=True)
    args = p.parse_args()
    return args

def _get_svtype_from_alt(alt):
    svtypes = {'DEL', 'DUP', 'INV', 'TRA'}
    if alt[0] == '<' and alt[-1] == '>':
        svtype = alt[1:-1]
    else:
        svtype = 'TRA'
    assert svtype in svtypes, svtype
    return svtype

def make_consensus_bedpe(df, sources, goi, consensus_cutoff=2):
    sv_cols = ['#chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
            'name', 'score', 'strand1', 'strand2', 'type', 'gene1', 'gene2']
    bedpe = pd.DataFrame(columns=sv_cols)
    for _, row in df.iterrows():
        chrom1 = row['chrom']
        end1 = int(row['pos'])
        start1 = end1 - 1
        chrom2, end2, strands = re.search('CHR2=([^;]+);.*END=(\d+);.*STRANDS=([\+\-]+);*', row['info']).groups()
        assert len(strands) == 2, strands
        #strand1, strand2 = strands
        end2 = int(end2)
        start2 = end2 - 1
        fmt = row['format'].split(':')
        cnt = 0
        ngoi = 0
        svids = []
        svtype = _get_svtype_from_alt(row['alt'])
        gene1s, gene2s = set(), set()
        for source in sources:
            gts = dict(zip(fmt, row[source].split(':')))
            if gts['ID'] != 'NaN':
                cnt += 1
                svid, gene1, gene2, strand1, strand2 = gts['ID'].split('__')
                gene1s.add(gene1); gene2s.add(gene2)
                ngoi += (gene1 in goi or gene2 in goi)
                svids.append(source)
        if ngoi > 0 or cnt >= consensus_cutoff:
            gene1, gene2 = '__'.join(list(gene1s)), '__'.join(list(gene2s))
            field = [chrom1, start1, end1, chrom2, start2, end2, '__'.join(svids), cnt, strand1, strand2, svtype, gene1, gene2]
            bedpe.loc[bedpe.shape[0]] = field
    return bedpe

if __name__ == "__main__":
    args = parse_args()
    vcf_cols = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format']
    vcf_cols += args.sources
    df = pd.read_table(args.vcf, names=vcf_cols, comment='#')

    consensus_cutoff = 2
    goi = {'MYC', 'BCL2', 'BCL6'}
    bedpe = make_consensus_bedpe(df, args.sources, goi=goi, consensus_cutoff=consensus_cutoff)
    bedpe.to_csv(args.bedpe, sep='\t', index=False)
