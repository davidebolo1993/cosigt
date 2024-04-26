rule bwamem2_index:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	'''
	input:
		rules.pyfaidx_extract.output
	output:
		multiext(config['output'] + '/pyfaidx/{region}.fasta', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')	
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.bwamem2_index.benchmark.txt'
	shell:
		'''
		bwa-mem2 index {input}
		'''

rule bwamem2_mem_samtools_sort:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	https://github.com/samtools/samtools
	'''
	input:
		ref=rules.pyfaidx_extract.output,
		idx=rules.bwamem2_index.output,
		sample=rules.samtools_fasta.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{region}.realigned.bam'
	threads:
		config['bwa-mem2']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{sample}.{region}.bwamem2_mem_samtools_sort.benchmark.txt'
	shell:
		'''
		bwa-mem2 mem \
		-t {threads} \
		{input.ref} \
		{input.sample} | samtools sort \
		-@ {threads} \
		- > {output}
		'''
