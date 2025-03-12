rule impg_project_batches:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.wfmash_align_batches.output,
		bed=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		config['output'] + '/impg/batches/{batch}.{region}.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/impg:0.2.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{batch}.{region}.impg_project.benchmark.txt'
	params:
		exclude='resources/extra/blacklist.txt'
	shell:
		'''
		impg \
		query \
		-p {input.paf} \
		-b {input.bed} | \
		grep -v \
		-E \
		-f {params.exclude} > {output}
		'''

rule concatenate_batches_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		expand(config['output'] + '/impg/batches/{{batch}}.{region}.bedpe', batch=sorted(batch_set))
	output:
		config['output'] + '/impg/{region}.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{region}.concatenate_batches_per_region.benchmark.txt'
	shell:
		'''
		for f in {input}; do
			cat $f >> {output}
		done
		'''