rule bedtools_merge:
	'''
	https://github.com/arq5x/bedtools2
	'''
	input:
		rules.concatenate_batches_per_region.output
	output:
		config['output'] + '/bedtools/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{region}.bedtools_merge.benchmark.txt'
	shell:
		'''
		bedtools sort \
		-i {input} | \
		bedtools merge \
		-d 200000 \
		-i - > {output}
		'''

rule filter_outliers:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.bedtools_merge.output
	output:
		config['output'] + '/bedtools/filtered/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{region}.filter_outliers.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/outliers.r \
		{input} \
		{output}
		'''