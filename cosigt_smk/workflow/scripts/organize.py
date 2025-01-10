#!/usr/bin/python3 env
#standard libraries
import os
import glob
import yaml
import argparse
from argparse import HelpFormatter

class CustomFormat(HelpFormatter):

	'''
	custom help format
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


def default_parameters(args):

	d=dict()
	#bwa-mem2
	d['bwa-mem2']=dict()
	d['bwa-mem2']['threads'] = args.aln_threads
	d['bwa-mem2']['mem_mb'] = args.aln_memory
	d['bwa-mem2']['time'] = args.aln_time
	#bwa-mem
	d['bwa']=dict()
	d['bwa']['threads'] = args.aln_threads
	d['bwa']['mem_mb'] = args.aln_memory
	d['bwa']['time'] =  args.aln_time
	#minimap2
	d['minimap2']=dict()
	d['minimap2']['threads'] = args.aln_threads
	d['minimap2']['mem_mb'] = args.aln_memory
	d['minimap2']['time'] =  args.aln_time
	d['minimap2']['preset'] =  args.aln_preset
	#samtools
	d['samtools']=dict()
	d['samtools']['threads'] = args.sam_threads
	d['samtools']['mem_mb'] = args.sam_memory
	d['samtools']['time'] =   args.sam_time
	#pggb
	d['pggb']=dict()
	d['pggb']['threads'] = args.pggb_threads
	d['pggb']['mem_mb'] = args.pggb_memory
	d['pggb']['time'] =  args.pggb_time
	d['pggb']['tmpdir'] = args.pggb_tmpdir
	d['pggb']['params'] =  args.pggb_params
	#wfmash
	d['wfmash']=dict()
	d['wfmash']['threads'] = args.wfmash_threads
	d['wfmash']['mem_mb'] = args.wfmash_memory
	d['wfmash']['time'] =  args.wfmash_time
	d['wfmash']['tmpdir'] = args.wfmash_tmpdir
	d['wfmash']['params'] =  args.wfmash_params
	#default
	d['default']=dict()
	d['default']['mem_mb'] = args.std_memory
	d['default']['time'] =  args.std_time
	#output
	d['output'] = args.output
	return d

def main():

	'''
	parse arguments and organize inputs for running cosigt properly without additional efforts from the users
	'''

	parser = argparse.ArgumentParser(prog='organize.py', description='''COsine SImilarity-based GenoTyper''', epilog='''Developed by Davide Bolognini @ Human Technopole''', formatter_class=CustomFormat) 
	#required
	required = parser.add_argument_group('Required I/O arguments')
	required.add_argument('-a', '--alignments', help='folder with read-level alignment files (.bam,.cram) - and their indexes (.bai/.csi,.crai). Individuals to genotype', metavar='FOLDER', required=True)
	required.add_argument('-r','--reference', help='reference file, .fasta format. Same reference the individuals to genotype are aligned to. Optionally, bgzip-compressed', metavar='FASTA', required=True)
	required.add_argument('--assemblies', help='assemblies, .fasta format. Optionally, bgzip-compressed', metavar='FASTA', required=True)
	required.add_argument('--roi', help='one or more regions of interest in .bed format, with chr in PanSN format (build#chr - e.g. grch38#chr1)', metavar='BED', required=True)
	#addition
	additional = parser.add_argument_group('Additional I/O arguments')
	additional.add_argument('--blacklist', help='assemblies (one per line) that should not be included in the analysis [None]', metavar='', required=False, default=None)
	additional.add_argument('--binds', help='additional paths to bind for singularity in /path/1,/path/2 format [/localscratch]', type=str, default='/localscratch')
	additional.add_argument('--tmp', help='SINGULARITY TMPDIR [/tmp]', type=str, default='/tmp')	
	additional.add_argument('--output', help='output folder [results]', metavar='FOLDER', default='results')
	additional.add_argument('--profile', help='use profile. If None, do not use profile and run on the local machine [config/slurm]', metavar='FOLDER', default='config/slurm', type=str)
	additional.add_argument('--samplemap', help='tsv file mapping each bam/cram basename to a user-defined id. If None, infer from bam/cram basename [None]', metavar='TSV', type=str, default=None)	
	additional.add_argument('--annotations', help='gene structures, .gtf format. Optionally gzip-compressed [None]', metavar='GTF', type=str, default=None)
	#metrics
	metrics = parser.add_argument_group('Specify #threads, memory and time requirements, temp directories')
	metrics.add_argument('--std_time', help='max time (minutes) - default [1]',type=int, default=1)
	metrics.add_argument('--std_memory', help='memory (mb) - default [500]',type=int, default=500)
	metrics.add_argument('--threads', help='number of concurrent cores for snakemake if no profile is provided (ignored otherwise) [1]', type=int, default=1)
	#alignment
	metrics.add_argument('--aln_threads', help='#threads - aligner [5]',type=int, default=5)
	metrics.add_argument('--aln_time', help='max time (minutes) - aligner [2]',type=int, default=5)
	metrics.add_argument('--aln_memory', help='max memory (mb) - aligner [5000]',type=int, default=5000)
	metrics.add_argument('--aln_preset', help='preset for minimap2 [map-ont] - ignore if not using the long branch of cosigt', type=str, default='map-ont')
	#samtools
	metrics.add_argument('--sam_threads', help='#threads - samtools (view) [2]',type=int, default=2)
	metrics.add_argument('--sam_time', help='max time (minutes) - samtools (view) [5]',type=int, default=5)
	metrics.add_argument('--sam_memory', help='max memory (mb) - samtools (view) [5000]',type=int, default=5000)
	#pggb
	metrics.add_argument('--pggb_threads', help='#threads - pggb [24]',type=int, default=24)
	metrics.add_argument('--pggb_time', help='max time (minutes) - pggb [35]',type=int, default=35)
	metrics.add_argument('--pggb_memory', help='max memory (mb) - pggb [30000]',type=int, default=30000)
	metrics.add_argument('--pggb_params', help='additional parameters for pggb [-c 2]',type=str, default='-c 2')
	metrics.add_argument('--pggb_tmpdir', help='temporary directory - pggb [working directory]',type=str, default=os.getcwd())
	#wfmash
	metrics.add_argument('--wfmash_threads', help='#threads - wfmash [24]',type=int, default=24)
	metrics.add_argument('--wfmash_time', help='max time (minutes) - wfmash [35]',type=int, default=35)
	metrics.add_argument('--wfmash_memory', help='max memory (mb) - wfmash [30000]',type=int, default=30000)
	metrics.add_argument('--wfmash_params', help='additional parameters for wfmash [-s 10k -p 95]',type=str, default='-s 10k -p 95')
	metrics.add_argument('--wfmash_tmpdir', help='temporary directory - wfmash [working directory]',type=str, default=os.getcwd())
	#parse args
	args = parser.parse_args()
	args.profile=None if args.profile == 'None' else args.profile
	#wd
	wd=os.getcwd()
	#default parameters
	d=default_parameters(args)
	#create all the output paths
	#config
	out_config_path=os.path.join(wd,'config')
	os.makedirs(out_config_path,exist_ok=True)
	out_yaml_tmp=os.path.join(out_config_path, 'config.yaml.tmp')
	out_yaml=os.path.join(out_config_path, 'config.yaml')
	out_samples=os.path.join(out_config_path, 'samples.tsv')
	#resources
	out_resources=os.path.join(wd,'resources')
	out_aln=os.path.join(out_resources, 'alignments')
	os.makedirs(out_aln, exist_ok=True)
	out_fasta=os.path.join(out_resources, 'assemblies')
	os.makedirs(out_fasta, exist_ok=True)
	out_ref=os.path.join(out_resources, 'reference')
	os.makedirs(out_ref, exist_ok=True)
	out_extra=os.path.join(out_resources, 'extra')
	os.makedirs(out_extra,exist_ok=True)
	blcklst_out=os.path.join(out_extra, 'blacklist.txt')
	out_regions=os.path.join(out_resources, 'regions')
	os.makedirs(out_regions,exist_ok=True)
	out_annotations=os.path.join(out_resources, 'annotations')
	os.makedirs(out_annotations,exist_ok=True)
	#blacklist of samples to exclude
	blcklst=[]
	if args.blacklist is not None:
		with open(args.blacklist, 'r') as bad_samples_in, open(blcklst_out, 'w') as bad_samples_out:
			for line in bad_samples_in:
				blcklst.append(line.rstrip())
				bad_samples_out.write(line)
	else:
		open(blcklst_out, 'w').close()
	#symlink alignments
	alns=sorted([x for x in glob.glob(args.alignments + '/**/*am*', recursive=True) if os.path.isfile(x)])
	samplesmap=dict()
	if args.samplemap is not None:
		with open(args.samplemap, 'r') as samples_in:
			for line in samples_in:
				sid,sname=line.rstrip().split('\t')
				samplesmap[sid]=sname
	with open(out_samples, 'w') as samples_out:
		samples_out.write('sample_id\talignment\n')
		for aln in alns:
			bnaln=os.path.basename(aln)
			out_aln_file=os.path.join(out_aln, bnaln)
			try:
				os.symlink(os.path.abspath(aln), out_aln_file) #error out if this exists
			except:
				pass #do not symlink again if exists
			if aln.endswith('am'): #this is not an index, rather a true alignment
				if args.samplemap is None:
					sample_name='.'.join(bnaln.split('.')[:-1])
					samples_out.write(sample_name + '\t' + out_aln_file + '\n')
				else:
					samples_out.write(samplesmap[bnaln] + '\t' + out_aln_file + '\n')
	#add to config
	d['samples'] = out_samples
	#annotations
	out_annotations_file = ''
	if args.annotations is not None:
		out_annotations_file=os.path.join(out_annotations, os.path.basename(args.annotations))
		#bind if needed
		args.binds += ',' + os.path.dirname(os.path.abspath(args.annotations))
	d['annotations'] = out_annotations_file
	#symlink assemblies 
	out_assemblies_file=os.path.join(out_fasta, os.path.basename(args.assemblies))
	try:
		os.symlink(os.path.abspath(args.assemblies), out_assemblies_file)
	except:
		pass
	#add to config
	d['assemblies'] = out_assemblies_file
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
	d['path'] = ''
	with open(args.roi) as bed_in:
		for line in bed_in:
			l=line.rstrip().split('\t')
			region=l[0].replace('#','_') + '_' + l[1] + '_' + l[2]
			d['path'] = l[0] if d['path'] == '' else d['path']
			#put regions in the config fille
			d['region'].append(region)
			#also write in the dedicated space
			region_out=os.path.join(out_regions, region+'.bed')
			with open(region_out, 'w') as out_region:
				out_region.write(l[0] + '\t' + l[1] + '\t' + l[2]+'\n')
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
	#write command - singularity
	singpath=','.join(list(set([os.path.abspath(args.alignments),os.path.dirname(os.path.abspath(args.assemblies)),os.path.abspath(args.output),os.path.dirname(os.path.abspath(args.reference)),args.binds,os.path.abspath(args.pggb_tmpdir),os.path.abspath(args.wfmash_tmpdir)])))
	if args.profile is not None:
		command_singularity_out='SINGULARITY_TMPDIR=' + os.path.abspath(args.tmp) + ' snakemake --profile ' + args.profile + ' --singularity-args "-B '+ singpath + ' -e" cosigt'
		with open('snakemake.singularity.profile.run.sh', 'w') as out:
			out.write(command_singularity_out + '\n')
		#write command - conda
		command_conda_out='snakemake --profile ' + args.profile + ' cosigt'
		with open('snakemake.conda.profile.run.sh', 'w') as out:
			out.write(command_conda_out + '\n')
	else:
		command_singularity_out='SINGULARITY_TMPDIR=' + os.path.abspath(args.tmp) + ' snakemake --use-singularity --singularity-args "-B '+ singpath + ' -e" -j ' + str(args.threads) + ' cosigt'
		with open('snakemake.singularity.run.sh', 'w') as out:
			out.write(command_singularity_out + '\n')
		#write command - conda
		command_conda_out='snakemake --use-conda -j ' + str(args.threads) + ' cosigt'
		with open('snakemake.conda.run.sh', 'w') as out:
			out.write(command_conda_out + '\n')

if __name__ == '__main__':
	
	main()
