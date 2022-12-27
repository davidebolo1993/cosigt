rule gafpack_coverage:
	'''
	gafpack
	'''
	input:
		gfa=rules.odgi_chop.output,
		gaf=rules.gfa_inject.output
	output:
		"results/cosigt_results/{sample}/{sample}.x.gafpack.gz"
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
