rule impg_project:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.wfmash_align.output,
		bed=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		config['output'] + '/impg/{region}.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{region}.impg_project.benchmark.txt'
	params:
		exclude='resources/extra/blacklist.txt'
	shell:
		'''
		impg \
		query \
		-p {input.paf} \
		-b {input.bed} | \
		grep -v \
		-f {params.exclude} > {output}
		'''
