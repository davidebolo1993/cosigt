rule pyfaidx_extract:
	'''
	https://github.com/mdshw5/pyfaidx/
	'''
	input:
		fasta=rules.add_target_to_queries.output,
		bed=rules.bedtools_merge.output
	output:
		config['output'] + '/pyfaidx/{region}.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.pyfaidx_extract.benchmark.txt'
	shell:
		'''
		faidx \
		--bed {input.bed} \
		{input.fasta} > {output}
		'''

rule samtools_faidx_index:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.pyfaidx_extract.output
	output:
		config['output'] + '/pyfaidx/{region}.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.samtools_faidx_index.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''
