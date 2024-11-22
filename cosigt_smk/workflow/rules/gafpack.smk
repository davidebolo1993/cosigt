rule gafpack_coverage:
	'''
	https://github.com/pangenome/gafpack
	'''
	input:
		gfa=rules.odgi_view.output,
		gaf=rules.gfainject_inject.output
	output:
		config['output'] + '/gafpack/{sample}/{region}.gafpack.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	conda:
		'../envs/gafpack.yaml'
	benchmark:
		'benchmarks/{sample}.{region}.gafpack_coverage.benchmark.txt'
	shell:
		'''
		gafpack \
		-g {input.gfa} \
		-a {input.gaf} \
		--len-scale | gzip > {output}
		'''
