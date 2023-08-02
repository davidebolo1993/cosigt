rule gafpack_coverage:
	'''
	gafpack
	'''
	input:
		gfa=rules.odgi_view.output,
		gaf=rules.gfa_inject.output
	output:
		'results/gafpack/{sample}/{region}.gafpack.gz'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		gafpack \
		-g {input.gfa} \
		-a {input.gaf} | pigz > {output}
		'''
