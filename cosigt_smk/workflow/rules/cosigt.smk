rule cosigt_genotype:
	'''
	cosigt genotype
	'''
	input:
		zpath=rules.odgi_paths_matrix.output,
		xpack=rules.gafpack_coverage.output,
		cluster=rules.cluster.output
	output:
		geno=config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv',
		combo=config['output'] + '/cosigt/{sample}/{region}/sorted_combos.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	params:
		prefix=config['output'] + '/cosigt/{sample}/{region}'
	shell:
		'''
		cosigt -p {input.zpath} -g {input.xpack} -c {input.cluster} -b resources/extra/blacklist.txt -o {params.prefix}
		'''