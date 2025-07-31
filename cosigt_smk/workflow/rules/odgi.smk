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
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
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
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
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

rule filter_nodes:
	'''
	https://github.com/davidebolo1993/cosigt 
	- This keeps all the node at the moment but in principle useful to filter out certain nodes
	'''	
	input:
		graph_cov=rules.odgi_paths.output,
		node_length=rules.get_nodes_length.output
	output:
		mask=config['output'] + '/odgi/paths/{chr}/{region}/{region}.mask.tsv',
		shared=temp(config['output'] + '/odgi/paths/{chr}/{region}/{region}.shared.tsv')
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.filter_nodes.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/filter.r \
		{input.graph_cov} \
		{input.node_length} \
		no_filter \
		{output.mask}
		'''	

rule odgi_dissimilarity:
	'''
	https://github.com/pangenome/odgi
	- Compute the dissimilarity of each path with respect to each other path
	'''
	input:
		og=rules.pggb_construct.output,
		mask=rules.filter_nodes.output.mask
	output:
		config['output'] + '/odgi/dissimilarity/{chr}/{region}/{region}.tsv.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
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
		-i {input.og} \
		-m {input.mask} \
		--all \
		--distances | gzip > {output}
		'''

rule make_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	- DBSCAN clustering based on dissimilarities between paths
	'''
	input:
		matrix=rules.odgi_dissimilarity.output,
		shared=rules.filter_nodes.output.shared
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
			{input.matrix} \
			{output} \
			automatic \
			$(cut -f 3 {input.shared} | tail -1) \
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
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
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