rule cosigt_genotype:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		zpath=rules.odgi_paths_matrix.output,
		xpack=rules.gafpack_coverage.output,
		cluster=rules.make_clusters.output
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
	benchmark:
		'benchmarks/{sample}.{region}.cosigt_genotype.benchmark.txt'
	shell:
		'''
		cosigt \
		-p {input.zpath} \
		-g {input.xpack} \
		-c {input.cluster} \
		-b resources/extra/blacklist.txt \
		-o {params.prefix}
		'''

rule cosigt_genotype2:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		zpath=rules.odgi_paths_matrix.output,
		xpack=rules.gafpack_coverage.output,
		cluster=rules.make_clusters2.output
	output:
		geno=config['output'] + '/cosigt2/{sample}/{region}/cosigt_genotype.tsv',
		combo=config['output'] + '/cosigt2/{sample}/{region}/sorted_combos.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	params:
		prefix=config['output'] + '/cosigt2/{sample}/{region}'
	benchmark:
		'benchmarks/{sample}.{region}.cosigt_genotype2.benchmark.txt'
	shell:
		'''
		cosigt \
		-p {input.zpath} \
		-g {input.xpack} \
		-c {input.cluster} \
		-b resources/extra/blacklist.txt \
		-o {params.prefix}
		'''