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

	# Set -s, the segment length, which must not exceed the shortest sequence.
	# Values are unchanged from the original expressions; they are just produced
	# as integers now. The original built them through float arithmetic and
	# emitted '-s 500.0', which pggb parses as a segment length of 0, and
	# '-s 0.0' for anything under 100 bp.
	if min_len <= 1000:
		# round down to the nearest 100, but never below the sequence itself
		s_val = (min_len // 100) * 100
		if s_val < 100:
			s_val = min_len
		if s_val > 0:
			flags.append(f'-s {s_val}')
	elif min_len <= 5000:
		# The original min(1000, floor(min_len/1000)*1000) is always 1000 here,
		# since floor(...) >= 1 in this branch. Kept as-is: raising it would
		# change the graphs. Revisit if a larger segment is wanted for 1-5 kb.
		flags.append('-s 1000')

	# Set -x
	if num_seq > 50:
		flags.append('-x auto')

	return " ".join(flags)



if ALLELE_SOURCE != 'graph':

	rule pggb_construct:
		'''
		https://github.com/pangenome/pggb
		- Build pangenome graph using pggb
		'''
		input:
			fasta=rules.bedtools_getfasta.output.fasta,
			fai=rules.bedtools_getfasta.output.fai
		output:
			graph=outpath("pggb/{chr}/{region}/{region}.og")
		threads:
			config['pggb']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['pggb']['runtime']
		container:
			'docker://ghcr.io/pangenome/pggb:20250423145743e25486'
		conda:
			'../envs/pggb.yaml'
		benchmark:
			'benchmarks/{chr}.{region}.pggb_construct.benchmark.txt'
		params:
			prefix=outpath("pggb/{chr}/{region}"),
			flags=lambda wildcards, input: config['pggb']['params'] + ' ' + parse_fai_and_set_flags(input.fai),
			tmpdir=tmp_path(config['pggb']['tmpdir'], "{chr}", "{region}", "{region}"),
			pansn=config['pansn_prefix']
		shell:
			'''
			if [ -d {params.prefix} ]; then
				rm -rf {params.prefix}
			fi
			mkdir -p {params.tmpdir}
			pggb \
				-i {input.fasta} \
				-o {params.prefix} \
				-t {threads} \
				-D {params.tmpdir} \
				-n $(wc -l < {input.fai}) \
				{params.flags} \
			&& odgi paths -i {params.prefix}/*smooth.final.og -L | grep {params.pansn} > {params.prefix}/ref_path.txt \
			&& odgi sort -i {params.prefix}/*smooth.final.og -Y -H {params.prefix}/ref_path.txt -o {output.graph} -C {params.tmpdir} \
			&& rm {params.prefix}/ref_path.txt
			rm -rf {params.tmpdir}
			'''

else:

	rule pggb_construct:
		'''
		https://github.com/pangenome/odgi
		- allele_source: graph -- stage a user-provided graph where pggb would
		  have written its own, so every rule downstream is unchanged
		- Accepts .og or plain GFA: odgi converts GFA on read, and the staged
		  copy is always .og, so nothing downstream sees the difference
		- Re-sorted by the reference path exactly as the pggb branch does with a
		  graph it builds, so both sources hand downstream the same layout
		- Its paths are extracted as the allele FASTA, written as stored and not
		  flipped: gfainject places realigned reads back onto this graph by path
		  name, so the FASTA has to match the graph base for base
		- Fails here, on the first step for the region, if the graph has no path
		  carrying the PanSN reference prefix
		'''
		input:
			custom_graph_path
		output:
			graph=outpath("pggb/{chr}/{region}/{region}.og"),
			fasta=outpath("pggb/{chr}/{region}/{region}.paths.fa")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
		container:
			'docker://ghcr.io/pangenome/pggb:20250423145743e25486'
		conda:
			'../envs/pggb.yaml'
		benchmark:
			'benchmarks/{chr}.{region}.pggb_construct.benchmark.txt'
		params:
			prefix=outpath("pggb/{chr}/{region}"),
			tmpdir=tmp_path(config['pggb']['tmpdir'], "{chr}", "{region}", "{region}"),
			pansn=config['pansn_prefix']
		shell:
			'''
			mkdir -p {params.prefix} {params.tmpdir}
			# Listed first and separately so that a graph odgi cannot read fails on
			# odgi's own error, instead of being misreported below as a graph that
			# merely lacks a reference path.
			odgi paths -i {input} -L > {params.prefix}/all_paths.txt
			grep {params.pansn} {params.prefix}/all_paths.txt > {params.prefix}/ref_path.txt || true
			if [ ! -s {params.prefix}/ref_path.txt ]; then
				echo "Graph {input} has no path with the reference prefix {params.pansn}." >&2
				echo "cosigt needs the reference path in every graph, named with pansn_prefix." >&2
				exit 1
			fi
			odgi sort -i {input} -Y -H {params.prefix}/ref_path.txt -o {output.graph} -C {params.tmpdir}
			odgi paths -i {output.graph} -f > {output.fasta}
			rm -f {params.prefix}/ref_path.txt {params.prefix}/all_paths.txt
			rm -rf {params.tmpdir}
			'''
