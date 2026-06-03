#!/usr/bin/env python3
"""
Generate a per-region VCF from cosigt genotype TSV files (one per sample)
and the FASTA FAI that lists all haplotypes for that region.
"""

import argparse
import csv
import os
import re
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


def numbered_haplotype_columns(fieldnames):
    columns = []
    for name in fieldnames:
        match = re.fullmatch(r"haplotype\.(\d+)", name)
        if match:
            columns.append((int(match.group(1)), name))
    return [name for _, name in sorted(columns)]


def read_last_genotype(tsv_path):
    """
    Return (sample_id, haplotype_names) from the last non-empty data row.
    cosigt TSV columns include haplotype.1 ... haplotype.N, where N is ploidy.
    """
    last = None
    with open(tsv_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames is None:
            sys.exit(f'ERROR: {tsv_path} is empty')
        hap_cols = numbered_haplotype_columns(reader.fieldnames)
        if not hap_cols:
            sys.exit(f'ERROR: no haplotype.N columns found in {tsv_path}')
        sample_col = '#sample.id'
        if sample_col not in reader.fieldnames:
            sample_col = 'sample' if 'sample' in reader.fieldnames else reader.fieldnames[0]
        for row in reader:
            if row and any(value not in (None, '') for value in row.values()):
                last = row
    if last is None:
        sys.exit(f'ERROR: {tsv_path} is empty')
    haplotypes = [last[column] for column in hap_cols]
    if any(hap in (None, '') for hap in haplotypes):
        sys.exit(f'ERROR: missing haplotype value in {tsv_path}')
    return last[sample_col], haplotypes


def main():
    args = parse_args()

    alleles, allele_index = parse_fai(args.fai, args.pansn)
    chrom, pos, end       = ref_coords(alleles[0])

    # ── collect genotypes ─────────────────────────────────────────────────
    samples   = []
    genotypes = {}
    for tsv in sorted(args.tsv):
        sample_id, haplotypes = read_last_genotype(tsv)
        samples.append(sample_id)
        gt = []
        for hap in haplotypes:
            idx = allele_index.get(hap)
            if idx is None:
                print(f'WARNING: {hap} not in FAI for {sample_id}; setting to .', file=sys.stderr)
                idx = '.'
            gt.append(str(idx))
        genotypes[sample_id] = '|'.join(gt)

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
                  'Description="Phased genotype as allele indices">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                  + '\t'.join(samples) + '\n')

        gt_cols = '\t'.join(genotypes.get(s, '.') for s in samples)
        out.write(f'{chrom}\t{pos}\t.\tN\t<HAPLOTYPE>\t.\t.\t{info}\tGT\t{gt_cols}\n')


if __name__ == '__main__':
    main()
