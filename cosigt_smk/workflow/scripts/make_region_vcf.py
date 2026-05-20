#!/usr/bin/env python3
"""
Generate a per-region VCF from cosigt genotype TSV files (one per sample)
and the FASTA FAI that lists all haplotypes for that region.
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--fai',     required=True,           help='Region FASTA .fai (haplotypes)')
    p.add_argument('--tsv',     required=True, nargs='+',help='cosigt TSV files (one per sample)')
    p.add_argument('--output',  required=True,           help='Output VCF path')
    p.add_argument('--pansn',   required=True,           help='PanSN reference prefix, e.g. grch38#1#chr6')
    return p.parse_args()

def parse_fai(fai_path, pansn_prefix):
    """
    Return (alleles, allele_index).
    alleles[0]  = reference path -> index 0
    alleles[1:] = all others, sorted alphabetically -> indices 1..N-1
    """
    names = []
    with open(fai_path) as fh:
        for line in fh:
            name = line.split('\t')[0]
            if name:
                names.append(name)

    ref_names   = [n for n in names if n.startswith(pansn_prefix)]
    other_names = sorted(n for n in names if not n.startswith(pansn_prefix))

    if not ref_names:
        sys.exit(f'ERROR: no entry in {fai_path} starts with "{pansn_prefix}"')

    alleles      = ref_names + other_names
    allele_index = {name: i for i, name in enumerate(alleles)}
    return alleles, allele_index


def ref_coords(ref_name):
    """
    'grch38#1#chr6:32484347-32603538' -> ('chr6', 32484347, 32603538)
    """
    chrom_coord = ref_name.split('#')[-1]
    chrom, span = chrom_coord.split(':')
    start, end  = span.split('-')
    return chrom, int(start), int(end)


def read_last_genotype(tsv_path):
    """
    Return (sample_id, hap1_name, hap2_name) from the last non-empty line.
    cosigt TSV columns: sample  hap1  hap2  haplogroup1  haplogroup2  score
    """
    last = None
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line:
                last = line
    if last is None:
        sys.exit(f'ERROR: {tsv_path} is empty')
    parts = last.split('\t')
    if len(parts) < 3:
        sys.exit(f'ERROR: unexpected format in {tsv_path}: {last!r}')
    return parts[0], parts[1], parts[2]


def main():
    args = parse_args()

    alleles, allele_index = parse_fai(args.fai, args.pansn)
    chrom, pos, end       = ref_coords(alleles[0])

    # ── collect genotypes ─────────────────────────────────────────────────
    samples   = []
    genotypes = {}
    for tsv in sorted(args.tsv):
        sample_id, hap1, hap2 = read_last_genotype(tsv)
        samples.append(sample_id)
        idx1 = allele_index.get(hap1)
        idx2 = allele_index.get(hap2)
        if idx1 is None:
            print(f'WARNING: {hap1} not in FAI for {sample_id}; setting to .', file=sys.stderr)
            idx1 = '.'
        if idx2 is None:
            print(f'WARNING: {hap2} not in FAI for {sample_id}; setting to .', file=sys.stderr)
            idx2 = '.'
        genotypes[sample_id] = f'{idx1}|{idx2}'

    # ── INFO ──────────────────────────────────────────────────────────────
    allele_map = ','.join(f'{i}={name}' for i, name in enumerate(alleles))
    info       = f'END={end};ALLELES={allele_map}'
    contig_line = f'##contig=<ID={chrom}>\n'

    # ── write VCF ─────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write(contig_line)
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,'
                  'Description="End coordinate of the genotyped region">\n')
        out.write('##INFO=<ID=ALLELES,Number=1,Type=String,'
                  'Description="Allele index to haplotype name: 0=ref,1=alt1,...">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,'
                  'Description="Phased diploid genotype as allele indices">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                  + '\t'.join(samples) + '\n')

        gt_cols = '\t'.join(genotypes.get(s, '.|.') for s in samples)
        out.write(f'{chrom}\t{pos}\t.\tN\t<HAPLOTYPE>\t.\t.\t{info}\tGT\t{gt_cols}\n')


if __name__ == '__main__':
    main()