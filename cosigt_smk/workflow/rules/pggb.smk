import math
import os

def parse_fai_and_set_flags(fai_input):
	'''
	https://github.com/davidebolo1993/cosigt
	- Dinamically set flags for pggb
	'''
	#this stuff should be namedlist, but need to check
	if hasattr(fai_input, 'fai'):
		fai_path = fai_input.fai
	elif isinstance(fai_input, (list, tuple)):
		fai_path = fai_input[0]  # just in case
	else:
		fai_path = fai_input  # assume string or os.PathLike

	if not isinstance(fai_path, (str, os.PathLike)):
		raise TypeError(f'Expected str or PathLike, got {type(fai_path)}')

	with open(fai_path) as f:
		lengths = [int(line.split('\t')[1]) for line in f]

	num_seq = len(lengths)
	total_len = sum(lengths)
	average_len = total_len / num_seq if num_seq > 0 else 0
	min_len = min(lengths)

	flags = []

	# Set -s
	if min_len <= 1000:
		s_val=round(min_len/1000,1)*1000
		flags.append(f'-s {s_val}')
	elif min_len <= 5000:
		s_val = min(1000, math.floor(min_len / 1000) * 1000)
		flags.append(f'-s {s_val}')

	# Set -x
	if num_seq > 50:
		flags.append('-x auto')

	return " ".join(flags)



rule pggb_construct:
	'''
	https://github.com/pangenome/pggb
	- Build pangenome graph using pggb
	'''
	input:
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai
	output:
		config['output'] + '/pggb/{chr}/{region}/{region}.og'
	threads:
		config['pggb']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	container:
		'docker://ghcr.io/pangenome/pggb:20250423145743e25486'
	conda:
		'../envs/pggb.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pggb_construct.benchmark.txt'
	params:
		prefix=config['output'] + '/pggb/{chr}/{region}',
		flags=lambda wildcards, input: config['pggb']['params'] + ' ' + parse_fai_and_set_flags(input.fai),
		tmpdir = config['pggb']['tmpdir'] + '/{chr}/{region}/{region}',
		pansn=config['pansn_prefix']
	shell:
		'''
		mkdir -p {params.tmpdir}
		pggb \
			-i {input.fasta} \
			-o {params.prefix} \
			-t {threads} \
			-D {params.tmpdir} \
			-n $(wc -l {input.fai}) \
			{params.flags} \
		&& odgi paths -i {params.prefix}/*smooth.final.og -L | grep {params.pansn} > {params.prefix}/ref_path.txt \
		&& odgi sort -i {params.prefix}/*smooth.final.og -Y -H {params.prefix}/ref_path.txt -o {output} \
		&& rm {params.prefix}/ref_path.txt
		rm -rf {params.tmpdir}
		'''
