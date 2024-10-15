rule cosigt_genotype:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		graph_cov_map=rules.odgi_paths_matrix.output,
		sample_cov_map=rules.gafpack_coverage.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv',
		config['output'] + '/cosigt/{sample}/{region}/sorted_combos.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	conda:
		'../envs/cosigt.yaml'
	params:
		prefix=config['output'] + '/cosigt/{sample}/{region}',
		sample_id='{sample}'
	benchmark:
		'benchmarks/{sample}.{region}.cosigt_genotype.benchmark.txt'
	shell:
		'''
		cosigt \
		-p {input.graph_cov_map} \
		-g {input.sample_cov_map} \
		-c {input.json} \
		-o {params.prefix} \
		-i {params.sample_id}
		'''