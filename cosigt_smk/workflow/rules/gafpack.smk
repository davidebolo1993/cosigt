rule gafpack_coverage:
	'''
	gafpack
	'''
	input:
		gfa=rules.odgi_view.output,
		gaf=rules.gfa_inject.output
	output:
		config['output'] + '/gafpack/{sample}/{region}.gafpack.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		gafpack \
		-g {input.gfa} \
		-a {input.gaf} | pigz > {output}
		'''