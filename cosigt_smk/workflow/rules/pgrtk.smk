rule pgrtk_bundle_decomposition:
	'''
	https://github.com/GeneDx/pgr-tk
	'''
	input:
		rules.samtools_faidx_extract.output
	output:
		config['output'] + '/pgrtk/pgr-pbundle-decomp/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/pgrtk:0.5.1'
	conda:
		'../envs/pgrtk.yaml'
	benchmark:
		'benchmarks/{region}.pgrtk_bundle_decomposition.txt'
	params:
		prefix=config['output'] + '/pgrtk/pgr-pbundle-decomp/{region}'
	shell:
		'''
		pgr-pbundle-decomp \
		{input} \
		{params.prefix}
		'''

rule build_bundles_structures:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.pgrtk_bundle_decomposition.output
	output:
		bundles_struct=config['output'] + '/pgrtk/bundles/{region}.bstruct.tsv',
		bundles_length=config['output'] + '/pgrtk/bundles/{region}.blength.tsv',
		bundles_table=config['output'] + '/pgrtk/bundles/{region}.all.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{region}.build_bundles_structures.txt'
	params:
		prefix=config['output'] + '/pgrtk/bundles/{region}'
	shell:
		'''
		Rscript \
		workflow/scripts/bundlestruct.r \
		{input} \
		{output.bundles_struct} \
		{output.bundles_length} \
		{output.bundles_table}
		'''

rule compute_bundles_distance:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		bundles_struct=rules.build_bundles_structures.output.bundles_struct,
		bundles_length=rules.build_bundles_structures.output.bundles_length
	output:
		config['output'] + '/pgrtk/bundles/{region}.bstruct.dist.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/pythonenv:3.13.3'
	conda:
		'../envs/python.yaml'
	benchmark:
		'benchmarks/{region}.compute_bundles_distance.txt'
	shell:
		'''
		python workflow/scripts/bundlesdist.py \
		{input.bundles_struct} \
		{input.bundles_length} \
		{output}
		'''

rule make_bundles_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.compute_bundles_distance.output
	output:
		config['output'] + '/pgrtk/cluster/{region}.clusters.json'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{region}.make_bundles_clusters.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/clusterbundles.r \
			{input} \
			{output} \
			automatic 
		'''

rule plot_bundles_clusters:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		tsv=rules.build_bundles_structures.output.bundles_table,
		json=rules.make_bundles_clusters.output
	output:
		config['output'] + '/pgrtk/plot/{region}.bclustered.pdf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{region}.plot_bundles_clusters.benchmark.txt'
	shell:
		'''
		Rscript \
			workflow/scripts/plotbundlescluster.r \
			{input.tsv} \
			{input.json} \
			{output}
		'''