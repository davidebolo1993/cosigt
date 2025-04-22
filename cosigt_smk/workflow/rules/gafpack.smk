rule gafpack_coverage:
	'''
	https://github.com/pangenome/gafpack
	'''
	input:
		gfa=rules.odgi_view.output,
		gaf=rules.gfainject_inject.output
	output:
		config['output'] + '/gafpack/{sample}/{chr}/{region}.gafpack.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/gafpack:0.1.2'
	conda:
		'../envs/gafpack.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.gafpack_coverage.benchmark.txt'
	shell:
		'''
		gafpack \
		--gfa {input.gfa} \
		--gaf {input.gaf} \
		--len-scale \
		--weight-queries | gzip > {output}
		'''
