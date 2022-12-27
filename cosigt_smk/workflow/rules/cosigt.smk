rule cosigt_genotype:
	'''
	cosigt genotype
	'''
	input:
		zpath=rules.odgi_build.output,
		xpack=rules.gafpack_coverage.output
	output:
		"results/cosigt_results/{sample}/best_genotype.tsv"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix="results/cosigt_results/{sample}"
	shell:
		'''
		cosigt {input.zpath} {input.xpack} {params.prefix}
		'''
