rule bedtools_merge:
	'''
	https://github.com/arq5x/bedtools2
	'''
	input:
		rules.impg_project.output
	output:
		config['output'] + '/bedtools/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.bedtools_merge.benchmark.txt'
	params:
		ref_region='resources/regions/{region}.bed'
	shell:
		'''
		bedtools sort \
		-i {input} | \
		bedtools merge \
		-d 100000 \
		-i - > {output} \
		&& cat {params.ref_region} >> {output}
		'''