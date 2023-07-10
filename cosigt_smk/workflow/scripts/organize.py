#!/usr/bin/python3

import os
import sys
import glob
import yaml
import shutil

def default_parameters():

    d=dict()

    #bwa-mem2
    d['bwa-mem2']=dict()
    d['bwa-mem2']['threads'] = 10
    d['bwa-mem2']['mem_mb'] = 10000
    d['bwa-mem2']['time'] = "00:30:00"

    #samtools
    d['samtools']=dict()
    d['samtools']['threads'] = 5
    d['samtools']['mem_mb'] = 5000
    d['samtools']['time'] = "00:30:00"

    return d


def main(cram_folder,graph_gfa,reference_fa,region_str,blacklist):

    d=default_parameters()

    out_config_path='config'
    out_yaml=os.path.join(out_config_path, 'config.yaml')
    out_samples=os.path.join(out_config_path, 'samples.tsv')

    #symlink crams
    out_resources_cram='resources/cram'
    os.makedirs(out_resources_cram, exist_ok=True)
    crams=sorted(glob.glob(cram_folder+'/*am*'))

    with open(out_samples, 'w') as outfile:

        outfile.write('sample_id\tcram\n')

        for cram in crams:

            if any(x in cram for x in blacklist):

                continue

            out_resouces_cram_file=os.path.join(out_resources_cram, os.path.basename(cram))
            os.symlink(cram, out_resouces_cram_file)
            
            if cram.endswith('am'):

                sample_name=os.path.basename(cram).split('.')[0]
                outfile.write(sample_name + '\t' + out_resouces_cram_file + '\n')


    #symlink graph gfa
    out_resources_graph='resources/graph'
    os.makedirs(out_resources_graph, exist_ok=True)
    out_resources_graph_file=os.path.join(out_resources_graph, os.path.basename(graph_gfa))
    os.symlink(graph_gfa, out_resources_graph_file)
    d['graph'] = out_resources_graph_file

    #symlink reference fasta
    out_resources_ref='resources/ref'
    os.makedirs(out_resources_ref, exist_ok=True)
    out_resources_ref_file=os.path.join(out_resources_ref, os.path.basename(reference_fa))
    os.symlink(reference_fa, out_resources_ref_file)
    d['reference'] = out_resources_ref_file

    #region
    d['region'] = region_str

    #samples
    d['samples'] = out_samples

    #dump
    yml_out=open(out_yaml, 'w')
    yaml.dump(d,yml_out)
    yml_out.close()

    #write blacklist
    os.makedirs('resources/extra', exist_ok=True)
        
    with open('resources/extra/bad_samples.txt', 'w') as outblck:

        for b in blacklist:

            outblck.write(b +'\n')


    #singularity
    with open('singularity_bind_paths.csv', 'w') as singpath:

        singpath.write(','.join([cram_folder,os.path.dirname(graph_gfa),os.path.dirname(reference_fa),'/localscratch']) + '\n')
    

if __name__ == '__main__':

    cram_folder=os.path.abspath(sys.argv[1])
    graph_gfa=os.path.abspath(sys.argv[2])
    reference_fa=os.path.abspath(sys.argv[3])
    region_str=sys.argv[4]
    
    try:

        blacklist_file=os.path.abspath(sys.argv[5])
        fin = open(blacklist_file)
        blacklist=[line.rstrip() for line in fin]
        fin.close()

    except:

        blacklist=[]

    main(cram_folder,graph_gfa,reference_fa,region_str,blacklist)