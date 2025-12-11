rule gafpack_coverage:
	'''
	https://github.com/pangenome/gafpack
	- Calculate the read coverage for each node in the graph
	'''
	input:
		gfa=rules.odgi_utils.output.gfa,
		gaf=rules.gfainject_inject.output
	output:
		config['output'] + '/gafpack/{sample}/{chr}/{region}/{region}.gafpack.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['mid']['time']
	container:
		'docker://davidebolo1993/gafpack:0.1.3'
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
