#!/usr/bin/python3 env

#standard libraries

import os
import sys
import glob
import yaml
import shutil
import argparse
from argparse import HelpFormatter


class CustomFormat(HelpFormatter):

	'''
	custo help format
	'''

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()



def double_quote(time):

	return '"%s"' % time

def default_parameters(args):

	d=dict()

	#bwa-mem2
	d['bwa-mem2']=dict()
	d['bwa-mem2']['threads'] = args.aln_threads
	d['bwa-mem2']['mem_mb'] = args.aln_memory
	d['bwa-mem2']['time'] = double_quote(args.aln_time)

	#bwa-mem
	d['bwa']=dict()
	d['bwa']['threads'] = args.aln_threads
	d['bwa']['mem_mb'] = args.aln_memory
	d['bwa']['time'] = double_quote(args.aln_time)

	#minimap2
	d['minimap2']=dict()
	d['minimap2']['threads'] = args.aln_threads
	d['minimap2']['mem_mb'] = args.aln_memory
	d['minimap2']['time'] = double_quote(args.aln_time)

	#samtools
	d['samtools']=dict()
	d['samtools']['threads'] = args.sam_threads
	d['samtools']['mem_mb'] = args.sam_memory
	d['samtools']['time'] =  double_quote(args.sam_time)

	#odgi
	d['odgi']=dict()
	d['odgi']['threads'] = args.odgi_threads
	d['odgi']['mem_mb'] = args.odgi_memory
	d['odgi']['time'] =  double_quote(args.odgi_time)

	return d


def main():

	'''
	parse arguments and organize inputs for running cosigt properly without additional efforts from the users
	'''

	parser = argparse.ArgumentParser(prog='cosigt', description='''COsine SImilarity-based GenoTyper''', epilog='''This program was developed by Davide Bolognini at Human Technopole.''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-r','--reference', help='reference genome in FASTA format', metavar='FASTA', required=True)
	required.add_argument('-g', '--graph', help='pangenome graph in GFA format', metavar='GFA', required=True)
	required.add_argument('-a', '--alignment', help='folder with alignment files (bam,cram) - and their index - to use', metavar='FOLDER', required=True)
	required.add_argument('--roi', help='one or more regions of interest in BED format', metavar='BED', required=True)

	additional = parser.add_argument_group('Additional I/O arguments')

	additional.add_argument('--blacklist', help='blacklist of samples (one per line) that should not be included in the analysis [None]', metavar='', required=False, default=None)
	additional.add_argument('--path', help='path in the graph to use as a reference [grch38]',type=str, default="grch38")

	metrics = parser.add_argument_group('Specify #threads, memory and time requirements')

	#alignment
	metrics.add_argument('--aln_threads', help='threads - aligner [10]',type=int, default=10)
	metrics.add_argument('--aln_time', help='max time (hh:mm:ss) - aligner ["00:30:00"]',type=str, default="00:30:00")
	metrics.add_argument('--aln_memory', help='max memory (mb) - aligner[10000]',type=int, default=10000)

	#samtools extraction and sort
	metrics.add_argument('--sam_threads', help='threads - samtools (view/sort) commands [5]',type=int, default=5)
	metrics.add_argument('--sam_time', help='max time (hh:mm:ss) - samtools (view/sort) commands ["00:10:00"]',type=str, default="00:15:00")
	metrics.add_argument('--sam_memory', help='max memory (mb) - samtools (view/sort) commands [5000]',type=int, default=5000)

	#odgi
	metrics.add_argument('--odgi_threads', help='threads - odgi (build/extract) commands [10]',type=int, default=10)
	metrics.add_argument('--odgi_time', help='max time (hh:mm:ss) - odgi (build/extract) commands ["00:30:00"]',type=str, default="00:30:00")
	metrics.add_argument('--odgi_memory', help='max memory (mb) - odgi (build/extract) commands [30000]',type=int, default=30000)


	args = parser.parse_args()

	#wd
	wd=os.getcwd()

	#default parameters
	d=default_parameters(args)


	#create all the output paths

	#config
	out_config_path=os.path.join(wd,'config')
	out_yaml_tmp=os.path.join(out_config_path, 'config.yaml.tmp')
	out_yaml=os.path.join(out_config_path, 'config.yaml')
	out_samples=os.path.join(out_config_path, 'samples.tsv')

	#resources
	out_resources=os.path.join(wd,'resources')
	out_aln=os.path.join(out_resources, 'alignment')
	os.makedirs(out_aln, exist_ok=True)
	out_graph=os.path.join(out_resources, 'graph')
	os.makedirs(out_graph, exist_ok=True)
	out_ref=os.path.join(out_resources, 'reference')
	os.makedirs(out_ref, exist_ok=True)
	out_extra=os.path.join(out_resources, 'extra')
	os.makedirs(out_extra,exist_ok=True)
	blcklst_out=os.path.join(out_extra, 'bad_samples.txt')
	out_regions=os.path.join(out_resources, 'regions')
	os.makedirs(out_regions,exist_ok=True)

	#sing
	out_sing=os.path.join(wd, 'singularity_bind_paths.csv')

	#blacklist of samples to exclude

	blcklst=[]

	if args.blacklist is not None:

		with open(args.blacklist, 'r') as bad_samples_in, open(blcklst_out, 'w') as bad_samples_out:

			for line in bad_samples_in:
				
				blcklst.append(line.rstrip())
				bad_samples_out.write(line)
		
	else:

		with open(blcklst_out, 'w') as bad_samples_out:

			pass

	#symlink alignments
	alns=sorted(glob.glob(args.alignment+'/*am*'))

	with open(out_samples, 'w') as samples_out:

		samples_out.write('sample_id\talignment\n')

		for aln in alns:

			if any(x.lower() in aln.lower() for x in blcklst):

				continue #skip blacklisted

			bnaln=os.path.basename(aln)
			out_aln_file=os.path.join(out_aln, bnaln)

			try:

				os.symlink(os.path.abspath(aln), out_aln_file)
			
			except:

				pass #do not symlink again if exists

			if aln.endswith('am'): #this is not an index, rather a true alignment

				sample_name=bnaln.split('.')[0]
				samples_out.write(sample_name + '\t' + out_aln_file + '\n')

	#add to config
	d['samples'] = out_samples

	#symlink graph 
	out_graph_file=os.path.join(out_graph, os.path.basename(args.graph))

	try:

		os.symlink(os.path.abspath(args.graph), out_graph_file)
	
	except:

		pass

	#add to config
	d['graph'] = out_graph_file

	#symlink reference
	out_reference_file=os.path.join(out_ref, os.path.basename(args.reference))

	try:

		os.symlink(os.path.abspath(args.reference), out_reference_file)
	
	except:

		pass

	#add to config
	d['reference'] = out_reference_file

	#add to config
	d['region'] = list()

	with open(args.roi) as bed_in:

		for line in bed_in:

			l=line.rstrip().split('\t')
			region=l[0] + ':' + l[1] + '-' + l[2]

			#put regions in the config fille
			d['region'].append(region)
		      
			#also write in the dedicated space
			region_out=os.path.join(out_regions, region+'.bed')

			with open(region_out, 'w') as out_region:

				out_region.write(args.path +'#' + l[0] + '\t' + l[1] + '\t' + l[2] + '\n')

	
	#dump config
	yml_out=open(out_yaml_tmp, 'w')
	yaml.dump(d,yml_out)
	yml_out.close()

	#remove single quotes

	with open(out_yaml_tmp) as filein, open(out_yaml, 'w') as fileout:

		for line in filein:
	
			line=line.replace("'","")
			fileout.write(line)

	os.remove(out_yaml_tmp)

	#write singularity paths for those using singularity
	with open(out_sing, 'w') as singpath:

		singpath.write(','.join([os.path.abspath(args.alignment),os.path.dirname(os.path.abspath(args.graph)),os.path.dirname(os.path.abspath(args.reference)),'/localscratch']) + '\n')


if __name__ == '__main__':
	
	main()