rule make_reference_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	- Get the reference .bed based on the original .bed
	'''
	input:
		lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		config['output'] + '/samtools/reference_bed/{chr}/{region}/{region}.bed.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.make_reference_bed.benchmark.txt'
	params:
		chrom='{chr}'
	shell:
		'''
		grep -w {params.chrom} {input} | gzip > {output}
		'''

rule make_alignment_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	- Make alignment .bed - that is, .bed for extracting reads from the original .cram files
	- This can be == to the reference bed if no alt contig where provided, can be different otherwise
	'''
	input:
		ref_bed=rules.make_reference_bed.output,
		ori_bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		config['output'] + '/samtools/alignment_bed/{chr}/{region}/{region}.bed.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_alignment_bed.benchmark.txt'
	shell:
		'''
		bedtools intersect \
			-a {input.ref_bed} \
			-b <(bedtools sort -i {input.ori_bed}) \
			-nonamecheck \
			-u | gzip > {output}
		bedtools intersect \
			-a  <(bedtools sort -i {input.ori_bed}) \
			-b {input.ref_bed} \
			-nonamecheck \
			-v | gzip >> {output}
		'''
