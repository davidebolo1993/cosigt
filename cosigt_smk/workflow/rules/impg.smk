rule impg_project_batches:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.wfmash_align_batches.output,
		bed=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		config['output'] + '/impg/batches/{region}/{batch}.bedpe'
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
		'benchmarks/{region}.{batch}.impg_project_batches.benchmark.txt'
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
	input:
		expand(config['output'] + '/impg/batches/{region}/{batch}.bedpe', region='{region}', batch=sorted(batch_set))
	output:
		config['output'] + '/impg/merged/{region}.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{region}.concatenate_batches_per_region.txt'
	shell:
		'''
		for f in {input}; do
			cat $f >> {output}
		done
		'''