rule impg_project:
	'''
	https://github.com/ekg/impg
	'''
	input:
		paf=config['paf'],
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		config['output'] + '/impg/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.impg_project.benchmark.txt'
	params:
		exclude='resources/extra/blacklist.txt'
	shell:
		'''
		impg \
		-p {input.paf} \
		-b {input.region} \
		-x | \
		cut -f 1-3 | \
		grep -v \
		-f {params.exclude} > {output}
		'''