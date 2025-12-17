#!/usr/bin/env python3
import os
import glob
import yaml
import argparse
import sys
import gzip


def make_default_config(tmp) -> dict:
    '''
    Make default config
    '''
    config=dict()
    #RESOURCES FOR THE MAIN BRANCH
    #THESE ARE DEFAULT RESOURCES USERS MAY WANT TO ADJUST IN THE CONFIG
    #config/config.yaml
    #bwamem2_index and bwamem2_mem_samtools_sort
    config['bwa-mem2']=dict()
    config['bwa-mem2']['threads'] = 5
    config['bwa-mem2']['mem_mb'] = 10000
    config['bwa-mem2']['time'] = 4
    #minimap2 for all-vs-all alignment
    config['minimap2']=dict()
    config['minimap2']['ava']=dict()
    config['minimap2']['ava']['threads'] = 4
    config['minimap2']['ava']['mem_mb'] = 5000
    config['minimap2']['ava']['time'] =  2
    #samtools for samtools_fasta_mapped rule in samtools.smk
    config['samtools']=dict()
    config['samtools']['fasta_mapped'] = dict()
    config['samtools']['fasta_mapped']['threads'] = 2
    config['samtools']['fasta_mapped']['mem_mb'] = 2000
    config['samtools']['fasta_mapped']['time'] = 2
    #samtools for samtools_fasta_unmapped rule in samtools.smk
    #depending on the size of the .bam this can be a lot slower than the fasta_mapped counterpart,
    #so we increment the resources here
    config['samtools']['fasta_unmapped'] = dict()
    config['samtools']['fasta_unmapped']['threads'] = 8
    config['samtools']['fasta_unmapped']['mem_mb'] = 8000
    config['samtools']['fasta_unmapped']['time'] = 8      
    #pggb
    #this really needs to be adjusted based on region length and parameters
    #but based on some benchmarking we did
    #on >300 loci in a single run
    config['pggb']=dict()
    config['pggb']['threads'] = 24
    config['pggb']['mem_mb'] = 20000
    config['pggb']['time'] =  40
    config['pggb']['tmpdir'] = tmp
    config['pggb']['params'] =  '-c 2 -k 101'
    #minimap2 for all-vs-one alignment
    #this depends a lot on the number of contigs
    #5 cores, 40 Gb, 40 min max - there are cases where we hit 80G(?) but overall
    #this should complete in minutes and do not exceed 10-20G
    config['minimap2']['avo']=dict()
    config['minimap2']['avo']['threads'] = 6
    config['minimap2']['avo']['mem_mb'] = 40000
    config['minimap2']['avo']['time'] =  50
    #meryl for db construction of the reference genome
    config['meryl']=dict()
    config['meryl']['threads'] = 10
    config['meryl']['mem_mb'] = 25000
    config['meryl']['time'] =  20
    #kfilt
    config['kfilt']=dict()
    config['kfilt']['threads'] = 8
    config['kfilt']['mem_mb'] = 20000
    config['kfilt']['time'] =  20
    #tiny rules use these resources instead
    #default - small
    config['default'] = dict()
    config['default']['small']=dict()
    config['default']['small']['mem_mb'] = 500
    config['default']['small']['time'] =  2
    #default - mid
    config['default']['mid']=dict()
    config['default']['mid']['mem_mb'] = 2000
    config['default']['mid']['time'] =  5
    #default - high
    config['default']['high']=dict()
    config['default']['high']['mem_mb'] = 10000
    config['default']['high']['time'] =  10
    #viz
    config['wally_viz'] = False
    config['svbyeye_viz'] = False
    config['odgi_viz'] = True
    config['pangene_viz'] = True
    #done
    print(f'Config template prepared!')
    return config

def validate_assembly(asm_path) -> bool:
    '''
    Validate assembly file
    '''
    if not os.path.exists(asm_path):
        print(f'Assembly file: {asm_path} does not exist!')
        return False
    if not os.access(asm_path, os.R_OK):
        print(f'Assembly file: {asm_path} is not readable!')
        return False
    if not asm_path.endswith('.fasta') and not asm_path.endswith('.fa') and not asm_path.endswith('.fasta.gz') and not asm_path.endswith('.fa.gz'):
        print(f'Assembly file: {asm_path} is not in the expected format (.fa/.fasta/.fa.gz/.fasta.gz)!')
        return False
    if not os.path.exists(asm_path + '.fai'):
        print(f'Assembly file: {asm_path} is not indexed - expected .fai; samtools faidx {asm_path} and retry!')
        return False
    if asm_path.endswith('.fasta.gz') or asm_path.endswith('.fa.gz'):
        if not os.path.exists(asm_path + '.gzi'):
            print(f'Assembly file: {asm_path} is not indexed - expected .gzi; samtools faidx {asm_path} and retry!')
            return False
    #read contigs and check they follow the PanSN specification
    with open(asm_path + '.fai', 'r') as idx_in:
        for line in idx_in:
            ctg_id=line.split('\t')[0]
            if len(ctg_id.split('#')) != 3:
                print(f'Contig: {ctg_id} in assembly file: {asm_path} does not follow PanSN-spec!')
                return False
    return True
    
def read_assemblies_file(assemblies_file) -> dict():
    '''
    Read and validate the assemblies specified in the tsv.
    Return a dict mapping each chromosome to its assembly
    '''
    assemblies_out=dict()
    if not os.path.exists(assemblies_file):
        print(f'Table with assemblies: {assemblies_file} does not exist!')
        sys.exit(1)
    if not os.access(assemblies_file, os.R_OK):
        print(f'Table with assemblies: {assemblies_file} is not readable!')
        sys.exit(1)
    with open(assemblies_file) as asm_in:
        for line in asm_in:
            entries=line.rstrip().split('\t')
            if len(entries) != 2:
                print(f'Table with assemblies: {assemblies_file} does not match the expected input format!')
                sys.exit(1)
            chrom=entries[0]
            asm=entries[1]
            if chrom in assemblies_out:
                print(f'Table with assemblies: {assemblies_file} contains duplicate chromosome entries!')
                sys.exit(1)
            if validate_assembly(asm):
                assemblies_out[chrom] = os.path.abspath(asm)
            else:
                sys.exit(1)
    print(f'Loaded table with assemblies {assemblies_file}!')
    return assemblies_out


# NEW: read_alignments_map replaces read_alignment_map + read_alignments(folder)
def read_alignments_map(alignment_map_file) -> dict():
    '''
    Read and validate alignments from a tab-separated file:
    alignment_path <TAB> sample_id
    Return a dict mapping each sample_id to its alignment path
    '''
    aln_dict = dict()
    if not os.path.exists(alignment_map_file):
        print(f'Alignments table: {alignment_map_file} does not exist!')
        sys.exit(1)
    if not os.access(alignment_map_file, os.R_OK):
        print(f'Alignments table: {alignment_map_file} is not readable!')
        sys.exit(1)

    with open(alignment_map_file, 'r') as map_in:
        for line in map_in:
            line = line.rstrip()
            if not line:
                continue
            try:
                aln_path, aln_id = line.split('\t')
            except ValueError:
                print(f'Alignments table: {alignment_map_file} has an invalid line (expected 2 columns): {line}')
                sys.exit(1)

            aln_path = os.path.abspath(aln_path)
            if aln_id in aln_dict:
                print(f'Duplicate alignment id {aln_id} in alignments table {alignment_map_file}!')
                sys.exit(1)

            if not validate_alignment(aln_path):
                sys.exit(1)

            aln_dict[aln_id] = aln_path

    print(f'Loaded alignments from {alignment_map_file}!')
    return aln_dict


def validate_alignment(alignment) -> bool:
    '''
    Validate alignment file
    '''
    if not os.path.exists(alignment):
        print(f'Alignment file: {alignment} does not exist!')
        return False
    if not os.access(alignment, os.R_OK):
        print(f'Alignment file: {alignment} is not readable!')
        return False
    if alignment.endswith('.bam'):
        if not os.path.exists(alignment + '.bai') and not os.path.exists(alignment + '.csi'):
            print(f'Alignment file: {alignment} is not indexed - expected .bai/.csi; samtools index {alignment} and retry!')
            return False
    if alignment.endswith('.cram'):
        if not os.path.exists(alignment + '.crai'):
            print(f'Alignment file: {alignment} is not indexed - expected .crai; samtools index {alignment} and retry!')
            return False
    return True


# OLD read_alignments removed


def read_bed(bed_file, asm_dict) -> dict():
    '''
    Read bed file and organize regions
    '''
    bed_dict=dict()
    if not os.path.exists(bed_file):
        print(f'Bed file: {bed_file} does not exist!')
        sys.exit(1)
    if not os.access(bed_file, os.R_OK):
        print(f'Bed file: {bed_file} is not readable!')
        sys.exit(1)
    with open(bed_file, 'r') as bedin:
        for line in bedin:
            bed_entry=line.rstrip().split('\t')
            if len(bed_entry) < 3:
                print(f'Invalid entry in {bed_file}: less then 3 columns for entry {bed_entry}')
            chrom=bed_entry[0]
            start=bed_entry[1]
            end=bed_entry[2]
            if len(bed_entry) >= 4:
                annot=bed_entry[3]
            else:
                annot="unknown"
            if len(bed_entry) == 5:
                alt=bed_entry[4]
            else:
               alt=None  
            if chrom not in asm_dict:
                print(f'Provided chromosome: {chrom} in bed file {bed_file} is missing in the assemblies')
                sys.exit(1)
            if chrom not in bed_dict:
                bed_dict[chrom] = [(chrom,start,end,annot,alt)]
            else:
                bed_dict[chrom].append((chrom,start,end,annot,alt))
    print(f'Loaded bed file {bed_file}!')    
    return bed_dict
    
def read_genome(genome_file) -> dict():
    '''
    Read and validate reference genome
    '''
    ref_dict=dict()
    if not os.path.exists(genome_file):
        print(f'Genome file: {genome_file} does not exist!')
        sys.exit(1)
    if not os.access(genome_file, os.R_OK):
        print(f'Genome file: {genome_file} is not readable!')
        sys.exit(1)
    if not genome_file.endswith('.fasta') and not genome_file.endswith('.fa') and not genome_file.endswith('.fasta.gz') and not genome_file.endswith('.fa.gz') and not genome_file.endswith('.fna') and not genome_file.endswith('.fna.gz'):
        print(f'Genome file: {genome_file} is not in the expected format (.fa/.fasta/.fna/.fa.gz/.fasta.gz/.fna.gz)!')
        sys.exit(1)
    if not os.path.exists(genome_file + '.fai'):
        print(f'Genome file: {genome_file} is not indexed - expected .fai; samtools faidx {genome_file} and retry!')
        sys.exit(1)
    if genome_file.endswith('.fasta.gz') or genome_file.endswith('.fa.gz') or genome_file.endswith('.fna.gz'):
        if not os.path.exists(genome_file + '.gzi'):
            print(f'Genome file: {genome_file} is not indexed - expected .gzi; samtools faidx {genome_file} and retry!')
            return False
    ref_dict[os.path.basename(genome_file)]=os.path.abspath(genome_file)
    print(f'Loaded reference file {genome_file}!')    
    return ref_dict


def read_gtf(gtf_file) -> dict():
    '''
    Read and validate reference genome
    '''
    gtf_dict=dict()
    if gtf_file is None:
        gtf_dict['NA'] = 'NA'
    else:
        if not os.path.exists(gtf_file):
            print(f'Gtf file: {gtf_file} does not exist!')
            sys.exit(1)
        if not os.access(gtf_file, os.R_OK):
            print(f'Gtf file: {gtf_file} is not readable!')
            sys.exit(1)
        if not gtf_file.endswith('.gtf') and not gtf_file.endswith('.gtf.gz'):
            print(f'Gtf file: {gtf_file} is not in the expected format (.gtf/.gtf.gz)!')
            sys.exit(1)
        gtf_dict[os.path.basename(gtf_file)]=os.path.abspath(gtf_file)
        print(f'Loaded gtf file {gtf_file}!')    
    return gtf_dict


def read_proteins(proteins_file) -> dict():
    '''
    Read and validate proteins sequences
    '''
    proteins_dict=dict()
    if proteins_file is None:
        proteins_dict['NA'] = 'NA'
    else:
        if not os.path.exists(proteins_file):
            print(f'Proteins FASTA file: {proteins_file} does not exist!')
            sys.exit(1)
        if not os.access(proteins_file, os.R_OK):
            print(f'Proteins FASTA file: {proteins_file} is not readable!')
            sys.exit(1)
        if not proteins_file.endswith('.fasta') and not proteins_file.endswith('.fa') and not proteins_file.endswith('.fasta.gz') and not proteins_file.endswith('.fa.gz') and not proteins_file.endswith('.fna') and not proteins_file.endswith('.fna.gz'):
            print(f'Proteins FASTA file: {proteins_file} is not in the expected format (.fa/.fasta/.fna/.fa.gz/.fasta.gz/.fna.gz)!')
            sys.exit(1)
        if not os.path.exists(proteins_file + '.fai'):
            print(f'Proteins FASTA file: {proteins_file} is not indexed - expected .fai; samtools faidx {proteins_file} and retry!')
            sys.exit(1)
        if proteins_file.endswith('.fasta.gz') or proteins_file.endswith('.fa.gz') or proteins_file.endswith('.fna.gz'):
            if not os.path.exists(proteins_file + '.gzi'):
                print(f'Proteins FASTA file: {proteins_file} is not indexed - expected .gzi; samtools faidx {proteins_file} and retry!')
                return False
        proteins_dict[os.path.basename(proteins_file)]=os.path.abspath(proteins_file)
        print(f'Loaded proteins FASTA file {proteins_file}!')    
    return proteins_dict


def validate_output(output_folder, config_yaml) -> dict():
    '''
    Validate output folder
    '''
    out_folder=os.path.abspath(output_folder)
    if not os.access(os.path.dirname(out_folder), os.W_OK):
        print(f'Parent folder of: {output_folder} is not writable!')
        sys.exit(1)
    config_yaml['SINGULARITY_BIND'].add(os.path.dirname(out_folder))
    config_yaml['output'] = out_folder
    print(f'Checked output folder {out_folder}!')
    return config_yaml


def write_assemblies(asm_dict, bed_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write assemblies to resources/assemblies
    Add chromosomes to analyse to config - not all provided if those are not necessary
    '''
    asm_dir=os.path.join(RESOURCES, 'assemblies')
    os.makedirs(asm_dir, exist_ok=True)
    config_yaml['chromosomes'] = set()

    for k,v in asm_dict.items():
        if k in bed_dict:
            asm_folder=os.path.join(asm_dir, k)
            os.makedirs(asm_folder, exist_ok=True)
            asm_name=os.path.basename(v)
            os.symlink(v, os.path.join(asm_folder,asm_name))
            os.symlink(v + '.fai', os.path.join(asm_folder,asm_name + '.fai'))
            if asm_name.endswith('.gz'):
                os.symlink(v + '.gzi', os.path.join(asm_folder,asm_name + '.gzi'))
            config_yaml['chromosomes'].add(k)
            config_yaml['SINGULARITY_BIND'].add(os.path.dirname(v))

    print(f'Added assemblies to {asm_dir}!')
    return config_yaml


def write_alignments(aln_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write alignments to resources/alignments
    Add samples to analyse to config
    '''
    aln_dir=os.path.join(RESOURCES, 'alignments')
    os.makedirs(aln_dir, exist_ok=True)
    config_yaml['samples'] = set()
    for k,v in aln_dict.items():
        if v.endswith('.bam'):
            aln_id=k + '.bam'
            aln_idx='.bai' if os.path.exists(v + '.bai') else '.csi'
        else:
            aln_id=k + '.cram'
            aln_idx='.crai'
        os.symlink(v, os.path.join(aln_dir,aln_id))
        os.symlink(v + aln_idx, os.path.join(aln_dir,aln_id + aln_idx))
        config_yaml['samples'].add(k)
        config_yaml['SINGULARITY_BIND'].add(os.path.dirname(v))
    print(f'Added alignments to {aln_dir}!')
    return config_yaml



def reference_contigs(config_yaml) -> dict:
    '''
    Given the reference index, return a dict
    with contig, start end
    '''
    ctg_dict=dict()
    ref_idx=config_yaml['reference'] + '.fai'
    with open(ref_idx, 'r') as refidx:
        for line in refidx:
            fields=line.rstrip().split('\t')
            chrom,start,end=fields[0],'0',fields[1]
            ctg_dict[chrom] = (chrom,start,end)
    return ctg_dict


def write_regions(bed_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write regions to resources/regions
    Add regions to config
    '''
    reg_dir=os.path.join(RESOURCES, 'regions')
    os.makedirs(reg_dir, exist_ok=True)
    config_yaml['regions'] = set()
    config_yaml['all_regions'] = os.path.abspath(os.path.join(reg_dir, 'all_regions.tsv'))
    contigs=reference_contigs(config_yaml)
    with open(config_yaml['all_regions'], 'w') as b_a_out:
        for k,v in bed_dict.items():
            bed_dir=os.path.join(reg_dir, k)
            os.makedirs(bed_dir, exist_ok=True)
            for subr in v:
                region_out='_'.join(subr[:-2])
                b_a_out.write('\t'.join(subr[:-1]) + '\n')
                bed_out=os.path.join(bed_dir, region_out + '.bed')
                with open(bed_out, 'w') as b_out:
                    b_out.write('\t'.join(subr[:-1]) + '\n')
                    if subr[-1] is not None:    
                        alts=subr[-1].split(',')
                        for alt in alts:
                            last=alt.rfind(':')
                            chr_alt,rest_alt=alt[:last], alt[last+1:]
                            start_alt,end_alt=rest_alt.split('-')
                            b_out.write('\t'.join([chr_alt,start_alt,end_alt,chr_alt]) + '\n')
                config_yaml['regions'].add(region_out)
    print(f'Added regions to {reg_dir}!')
    return config_yaml


def write_reference(genome_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write reference to resources/reference
    Add reference to config
    '''
    ref_dir=os.path.join(RESOURCES, 'reference')
    os.makedirs(ref_dir, exist_ok=True)
    for k,v in genome_dict.items():
        os.symlink(v, os.path.join(ref_dir,k))
        os.symlink(v + '.fai', os.path.join(ref_dir,k + '.fai'))
        if v.endswith('gz'):
            os.symlink(v + '.gzi', os.path.join(ref_dir,k + '.gzi'))
    config_yaml['reference'] = os.path.join(os.path.abspath(ref_dir),k)
    config_yaml['SINGULARITY_BIND'].add(os.path.dirname(v))
    print(f'Added reference to {ref_dir}!')
    return config_yaml


def write_gtf(gtf_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write gtf to resources/annotations
    Add gtf to config
    '''
    gtf_dir=os.path.join(RESOURCES, 'annotations')
    os.makedirs(gtf_dir, exist_ok=True)
    for k,v in gtf_dict.items():
        if k == 'NA':
            config_yaml['gtf'] = k
            break
        else:
            os.symlink(v, os.path.join(gtf_dir,k))
            config_yaml['gtf'] = os.path.join(os.path.abspath(gtf_dir),k)
            config_yaml['SINGULARITY_BIND'].add(os.path.dirname(v))
    print(f'Added gtf to {gtf_dir}!')
    return config_yaml


def write_proteins(proteins_dict, config_yaml, RESOURCES) -> dict:
    '''
    Write proteins fasta to resources/annotations
    Add proteins to config   
    '''
    proteins_dir=os.path.join(RESOURCES, 'annotations')
    os.makedirs(proteins_dir, exist_ok=True)
    for k,v in proteins_dict.items():
        if k == 'NA':
            config_yaml['proteins'] = k
            break
        else:
            os.symlink(v, os.path.join(proteins_dir,k))
            os.symlink(v + '.fai', os.path.join(proteins_dir,k + '.fai'))
            if v.endswith('gz'):
                os.symlink(v + '.gzi', os.path.join(proteins_dir,k + '.gzi'))
            config_yaml['proteins'] = os.path.join(os.path.abspath(proteins_dir),k)
            config_yaml['SINGULARITY_BIND'].add(os.path.dirname(v))
    print(f'Added proteins to {proteins_dir}!')
    return config_yaml

def validate_flagger(flagger_blacklist, config_yaml, RESOURCES):
    '''
    Write flagger blacklist
    Add flagger blacklist to config
    '''
    flagger_blacklist_dir=os.path.join(RESOURCES, 'flagger')
    os.makedirs(flagger_blacklist_dir, exist_ok=True)
    flagger_blacklist_out=os.path.join(flagger_blacklist_dir, 'flagger_blacklist.bed')
    if flagger_blacklist is None:
        open(flagger_blacklist_out, 'w').close()
    else:
        if not os.path.exists(flagger_blacklist):
            print(f'Flagger blacklist file: {flagger_blacklist} does not exist!')
            sys.exit(1)
        if not os.access(flagger_blacklist, os.R_OK):
            print(f'Flagger blacklist file: {flagger_blacklist} is not readable!')
            sys.exit(1)
        with open(flagger_blacklist, 'r') as fb_in, open(flagger_blacklist_out, 'w') as fb_out:
            for line in fb_in:
                ctg_id=line.split('\t')[0]
                if len(ctg_id.split('#'))!=3:
                    print(f'Flagger blacklist: {flagger_blacklist} contains contig {ctg_id} which does not follow PanSN-spec!')
                    sys.exit(1)
                fb_out.write(line)
    print(f'Wrote flagger blacklist file to {flagger_blacklist_out}!')    
    config_yaml['flagger_blacklist'] = os.path.abspath(flagger_blacklist_out)
    return config_yaml

def write_config(config_yaml, config_out):
    '''
    Write config file
    '''
    yml_out=open(config_out, 'w')
    yaml.dump(config_yaml,yml_out)
    yml_out.close()
    print(f'Wrote config to {config_out}!')

def find_optimal_bindings(paths, min_coverage_threshold=2):
    """
    Find minimal number of paths to bind for singularity
    """
    if not paths:
        return []

    dir_coverage = {}
    for path in paths:
        components = path.split('/')
        current_path = ""
        for component in components:
            if not component:
                continue
            if current_path:
                current_path = f"{current_path}/{component}"
            else:
                current_path = f"/{component}"
            if current_path not in dir_coverage:
                dir_coverage[current_path] = 0
            dir_coverage[current_path] += 1
    path_bindings = {}
    for path in paths:
        components = path.split('/')
        best_binding = path        
        current_path = ""
        for i, component in enumerate(components):
            if not component:
                continue
            if current_path:
                current_path = f"{current_path}/{component}"
            else:
                current_path = f"/{component}"            
            if dir_coverage[current_path] >= min_coverage_threshold:
                best_binding = current_path
        path_bindings[path] = best_binding
    bindings = set(path_bindings.values())
    # Optimize the bindings - if a binding is fully contained within another, remove it
    optimized_bindings = set()
    for binding in sorted(bindings, key=len):
        # Check if any of the existing optimized bindings contain this binding
        if any(binding.startswith(other_binding + '/') for other_binding in optimized_bindings):
            continue
        optimized_bindings.add(binding)
    if '/' in optimized_bindings and len(optimized_bindings) > 1:
        optimized_bindings.remove('/')    
    return sorted(optimized_bindings, key=len)

def setup_arg_parser():
    '''
    Simple argument parser
    '''
    parser = argparse.ArgumentParser(
        prog='organize.py', 
        description='COsine SImilarity-based GenoTyper', 
        epilog='Developed by Davide Bolognini @ Human Technopole', 
    )
    # Required arguments
    required = parser.add_argument_group('Required I/O arguments')
    required.add_argument('-a', '--assemblies', help='assemblies individuals to -r will be genotyped against. This is a tab-separated file mapping chromosomes in -b to a FASTA with contigs for that chromosome. FASTA can be bgzip-compressed and must be indexed', metavar='FASTA', required=True)
    required.add_argument('-g', '--genome', help='reference genome. This is the FASTA regions to -b refers to. FASTA can be bgzip-compressed and must be indexed', metavar='FASTA', required=True)
    # CHANGED: --reads is now a TSV map, not a folder
    required.add_argument(
        '-r', '--reads',
        help='individuals to genotype. This is a tab-separated file mapping each alignment file (1st column, BAM/CRAM, indexed) to an id (2nd column)',
        metavar='TSV',
        required=True
    )
    required.add_argument('-b', '--bed', help='regions to genotype. A standard 3-column BED file, but can have a 4th column to label the region and a 5th column listing comma-separated alternative contigs for that region', metavar='BED', required=True)
    required.add_argument('-o', '--output', help='output folder. This is where results from cosigt pipeline will be stored', metavar='FOLDER', required=True)
    optional = parser.add_argument_group('Optional arguments')
    # REMOVED: --map
    optional.add_argument('--gtf', help='Gene annotation on the reference chromosomes in GTF format [None]', metavar='', required=False, default=None)
    optional.add_argument('--proteins', help='Protein-coding transcript translation sequences in FASTA format.  FASTA can be bgzip-compressed and must be indexed [None]', metavar='', required=False, default=None)
    optional.add_argument('--tmp', help='tmp directory. Will be used by singularity and pggb [/tmp]', metavar='', required=False, default='/tmp')
    optional.add_argument('--pansn', help='PanSN prefix naming for the reference genome [grch38#1#]', metavar='', required=False, default='grch38#1#')
    optional.add_argument('--profile', help='snakemake profile, if available [None]', metavar='', required=False, default=None)
    optional.add_argument('--flagger', help='regions to exclude for each contigs. This is a standard BED file coming from flagger, with contigs names matching those of the assemblies in -a [None]', metavar='', required=False, default=None)
    optional.add_argument('--conda', help='prepare for running using conda instead of singularity [False]', action='store_true')
    optional.add_argument('--threads', help='run snakemake using that many cores - ignored if using a profile [32]', metavar='', required=False, default=32)
    optional.add_argument('--no_pangene', help='DO NOT visualize coding genes on input haplotypes [False]', action='store_false')
    optional.add_argument('--no_odgi', help='DO NOT visualize node coverage on input haplotypes in a "odgi viz"-like format [False]', action='store_false')
    optional.add_argument('--wally', help='run wally to visualize reads-to-haplotypes realignment [False]', action='store_true')
    optional.add_argument('--svbyeye', help='run SVByEye to visualize predicted haplotypes-to-reference realignment [False]', action='store_true')

    return parser


def main():
    '''
    Maind    '''
    parser = setup_arg_parser()
    args = parser.parse_args()


    CWD=os.path.realpath(__file__)
    BASE=os.path.dirname(os.path.dirname(os.path.dirname(CWD)))
    RESOURCES=os.path.join(BASE, 'resources')
    os.makedirs(RESOURCES, exist_ok=True)
    CONFIG=os.path.join(BASE, 'config')
    os.makedirs(CONFIG, exist_ok=True)

    #check conflicts
    if (args.gtf is None and args.proteins is not None) or (args.gtf is not None and args.proteins is None):
        print(f"When providing argument to --gtf, also an argument to --proteins must be specified")
        sys.exit(1)


    if os.path.isdir(RESOURCES):
        if os.listdir(RESOURCES):
            print(f"Directory {RESOURCES} is not empty: clean it and retry!")
            sys.exit(1)


    #READ
    #initial config
    #config
    config=make_default_config(os.path.abspath(args.tmp))
    config['wally_viz'] = args.wally
    config['svbyeye_viz'] = args.svbyeye
    if args.svbyeye:
        print(f"When --svbyeye is specified, --conda must be set to False since SVByEye is not implemented in a dedicated conda environment")
        sys.exit(1)
    config['pangene_viz'] = args.no_pangene
    config['odgi_viz'] = args.no_odgi
    config['SINGULARITY_BIND'] = set()
    config['SINGULARITY_BIND'].add(os.path.abspath(args.tmp))
    config['pansn_prefix'] = args.pansn
    #assemblies
    asm_dict=read_assemblies_file(os.path.abspath(args.assemblies))
    # bed file
    bed_dict=read_bed(os.path.abspath(args.bed),asm_dict)
    #genome
    genome_dict=read_genome(os.path.abspath(args.genome))
    #gtf / proteins
    gtf_dict=read_gtf(args.gtf)
    proteins_dict=read_proteins(args.proteins)
    # NEW: read alignments from TSV map
    aln_dict=read_alignments_map(os.path.abspath(args.reads))
    #output
    config=validate_output(os.path.abspath(args.output), config)


    #WRITE    
    #assemblies
    config=write_assemblies(asm_dict, bed_dict, config, RESOURCES)
    #alignments
    config=write_alignments(aln_dict, config, RESOURCES)
    #reference
    config=write_reference(genome_dict, config, RESOURCES)    
    #regions
    config=write_regions(bed_dict, config, RESOURCES)
    #annotations
    config=write_gtf(gtf_dict, config, RESOURCES)
    config=write_proteins(proteins_dict, config, RESOURCES)
    #flagger blacklist
    config=validate_flagger(args.flagger, config, RESOURCES)
    
    #common paths to BIND for singularity
    paths=find_optimal_bindings(list(config['SINGULARITY_BIND']))
    del config['SINGULARITY_BIND']
    #write config
    config['chromosomes'] = list(config['chromosomes'])
    config['samples'] = list(config['samples'])
    config['regions'] = list(config['regions'])
    write_config(config, os.path.join(CONFIG, 'config.yaml'))

    #write cosigt command
    cmd='#!/bin/bash\n'
    if args.profile is not None:
        if not args.conda:
            cmd +='SINGULARITY_TMPDIR=' + os.path.abspath(args.tmp) + ' snakemake --profile ' + args.profile + ' --singularity-args "-B '+ ','.join(paths) + ' -e" cosigt'
        else:
            cmd += 'snakemake --profile ' + args.profile + ' cosigt'
    else: #no profile
        if not args.conda:
            cmd +='SINGULARITY_TMPDIR=' + os.path.abspath(args.tmp) + ' snakemake --use-singularity --singularity-args "-B '+ ','.join(paths) + ' -e" -j ' + str(args.threads) + ' cosigt'
        else:
            cmd += 'snakemake --use-conda -j ' + str(args.threads) + ' cosigt'

    cmd+=' --rerun-triggers=mtime --rerun-incomplete --scheduler greedy\n'
    cmd_out=os.path.join(BASE, 'cosigt_smk.sh')
    with open(cmd_out, 'w') as fout:
        fout.write(cmd)
    print(f'Wrote cosigt command to run the the pipeline in: {cmd_out}!')
    sys.exit(0)


if __name__ == "__main__":
    main()
