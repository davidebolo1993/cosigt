rule cosigt_genotype:
	'''
	cosigt genotype
	'''
	input:
		zpath=rules.odgi_paths_matrix.output,
		xpack=rules.gafpack_coverage.output
	output:
		geno='results/cosigt/{sample}/{region}/best_genotype.tsv',
		combo='results/cosigt/{sample}/{region}/combos.tsv'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix='results/cosigt/{sample}/{region}'
	shell:
		'''
		cosigt -p {input.zpath} -g {input.xpack} -b resources/extra/blacklist.txt -o {params.prefix}
		'''

