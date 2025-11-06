rule odgi_view:
	'''
	https://github.com/pangenome/odgi
	- Convert .og to .gfa - this is mainly necessary downstream
	'''
	input:
		rules.pggb_construct.output
	output:
		temp(config['output'] + '/odgi/view/{chr}/{region}/{region}.gfa.gz')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1753347183'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_view.benchmark.txt'
	shell:
		'''
		odgi view \
		-i {input} \
		-g | gzip > {output}
		'''

rule odgi_paths:
	'''
	https://github.com/pangenome/odgi
	- Compute paths coverage over each node in the graph
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/paths/{chr}/{region}/{region}.tsv.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1753347183'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_paths.benchmark.txt'
	shell:
		'''
		odgi paths \
		-i {input} \
		-H | \
		cut -f 1,4- | gzip > {output}
		'''

rule get_nodes_length:
	'''
	https://github.com/pangenome/odgi
	- Just compute the length of each node
	'''
	input:
		rules.odgi_view.output
	output:
		config['output'] + '/odgi/view/{chr}/{region}/{region}.node.length.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.get_node_length.benchmark.txt'
	shell:
		'''
		zgrep '^S' {input} | \
		awk '{{print("node."$2,length($3))}}' OFS="\\t" > {output}
		'''

rule panplexity_filter:
	'''
	https://github.com/davidebolo1993/cosigt 
	- This builds a mask for nodes, excluding those that are low-complexity
	'''	
	input:
		rules.odgi_view.output
	output:
		temp(config['output'] + '/panplexity/{chr}/{region}/{region}.mask.tsv')
	threads:
		4
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/panplexity:0.1.1'
	conda:
		'../envs/panplexity.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.panplexity_filter.benchmark.txt'
	shell:
		'''
		panplexity \
			--input-gfa {input} \
			-t auto \
			-k 16 \
			-w 100 \
			-d 100 \
			--complexity linguistic \
			-m {output} \
			--threads 4
		'''	

rule filter_nodes:
	'''
	https://github.com/davidebolo1993/cosigt
	- Additional mask for nodes, mask nodes exhibiting coverage spikes on the paths
	- Output a combined mask (panplexity + coverage) and some statistics/plots
	'''
	input:
		paths=rules.odgi_paths.output,
		mask=rules.panplexity_filter.output,
		lengths=rules.get_nodes_length.output
	output:
		config['output'] + '/odgi/paths/{chr}/{region}/{region}.mask.tsv'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.filter_nodes.benchmark.txt'
	params:
		prefix=config['output'] + '/odgi/paths/{chr}/{region}/{region}'
	shell:
		'''
		Rscript \
			workflow/scripts/coverage_outliers.r \
			{input.paths} \
			{params.prefix} \
			{input.lengths} \
			{input.mask}
		'''		

rule odgi_dissimilarity:
	'''
	https://github.com/pangenome/odgi
	- Compute the dissimilarity of each path with respect to each other path
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/dissimilarity/{chr}/{region}/{region}.tsv.gz'
	threads:
		5
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1753347183'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_dissimilarity.benchmark.txt'
	shell:
		'''
		odgi similarity \
		-i {input} \
		--all \
		--distances \
		-t {threads} | gzip > {output}
		'''

rule make_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	- DBSCAN clustering based on dissimilarities between paths
	'''
	input:
		rules.odgi_dissimilarity.output
	output:
		config['output'] + '/cluster/{chr}/{region}/{region}.clusters.json'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_clusters.benchmark.txt'
	params:
		threshold_file=config['output'] + '/odgi/paths/{chr}/{region}/{region}.shared.tsv'
	shell:
		'''
		Rscript workflow/scripts/cluster.r \
			{input} \
			{output} \
			automatic \
			100.0 \
			1
		'''

rule viz_odgi:
	'''
	https://github.com/davidebolo1993/cosigt
	- Visualize the node coverage using the clustering information
	- This serves as a (better) alternative to usual odgi viz
	'''
	input:
		graph_cov=rules.odgi_paths.output,
		json=rules.make_clusters.output,
		nodes_length=rules.get_nodes_length.output
	output:
		config['output'] + '/odgi/viz/{chr}/{region}/{region}.viz.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.viz_odgi.benchmark.txt'
	params:
		tsv=config['output'] + '/cluster/{chr}/{region}/{region}.clusters.medoids.tsv'
	shell:
		'''
		Rscript \
			workflow/scripts/viz_odgi.r \
			{input.graph_cov} \
			{input.nodes_length} \
			{input.json} \
			{output} \
			{params.tsv}
		'''
