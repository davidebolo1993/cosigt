rule odgi_chop:
	'''
	https://github.com/pangenome/odgi
	skipping this for the time being
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/chop/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_chop.benchmark.txt'
	shell:
		'''
		odgi chop  \
		-i {input} \
		-c 32 \
		-o {output}
		'''

rule odgi_view:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/view/{region}.gfa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_view.benchmark.txt'
	shell:
		'''
		odgi view \
		-i {input} \
		-g > {output}
		'''

rule odgi_paths_matrix:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/paths/matrix/{region}.tsv.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_paths_matrix.benchmark.txt'
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
	'''
	input:
		rules.odgi_view.output
	output:
		config['output'] + '/odgi/view/{region}.node.length.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{region}.get_node_length.benchmark.txt'
	shell:
		'''
		grep '^S' {input} | \
		awk '{{print("node."$2,length($3))}}' OFS="\\t" > {output}
		'''

rule filter_nodes:
	'''
	https://github.com/davidebolo1993/cosigt 
	'''	
	input:
		graph_cov=rules.odgi_paths_matrix.output,
		node_length=rules.get_nodes_length.output
	output:
		config['output'] + '/odgi/paths/matrix/{region}.mask.tsv'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{region}.filter_nodes.benchmark.txt'
	params:
		mask=config['mask']
	shell:
		'''
		Rscript workflow/scripts/filter.r \
		{input.graph_cov} \
		{input.node_length} \
		{params.mask} \
		{output}
		'''	

rule odgi_similarity:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		mask=rules.filter_nodes.output
	output:
		config['output'] + '/odgi/similarity/{region}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_similarity.benchmark.txt'
	shell:
		'''
		odgi similarity \
		-i {input.og} \
		-m {input.mask} > {output}
		'''

rule odgi_dissimilarity:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		mask=rules.filter_nodes.output
	output:
		config['output'] + '/odgi/dissimilarity/{region}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_dissimilarity.benchmark.txt'
	shell:
		'''
		odgi similarity \
		-i {input.og} \
		-m {input.mask} \
		--distances > {output}
		'''

rule make_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.odgi_similarity.output
	output:
		config['output'] + '/cluster/{region}.clusters.json'
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
		'benchmarks/{region}.make_clusters.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/cluster.r \
			{input} \
			{output}
		'''

rule make_dbscan_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.odgi_dissimilarity.output
	output:
		config['output'] + '/cluster_dbscan/{region}.clusters.json'
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
		'benchmarks/{region}.make_clusters.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/cluster_dbscan.r \
			{input} \
			{output} \
			0.95
		'''

rule odgi_viz:
	'''
	https://github.com/pangenome/odgi
	this is going to change once we can cluster odgi viz directly
	'''
	input:
		og=rules.pggb_construct.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/odgi/viz/{region}.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_viz.benchmark.txt'
	params:
		tsv=config['output'] + '/cluster/{region}.clusters.tsv'
	shell:
		'''
		odgi viz \
		-i {input.og} \
		-p <(cut -f 1 {params.tsv} | tail -n+2) \
		-m \
		-o {output}
		'''

rule subset_gtf:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		gtf=config['annotations'],
		bed=rules.make_reference_bed.output
	output:
		config['output'] + '/annotations/{region}.gtf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{region}.subset_gtf.benchmark.txt'
	shell:
		'''
		bedtools intersect \
		-wa \
		-a {input.gtf} \
		-b {input.bed} > {output}
		'''

rule make_annotation_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.subset_gtf.output
	output:
		config['output'] + '/annotations/{region}.bed'
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
		'benchmarks/{region}.make_annotation_bed.benchmark.txt'
	params:
		path=config['path']
	shell:
		'''
		Rscript workflow/scripts/annotate.r \
		{input} \
		{params.path} \
		{output}
		'''		

rule odgi_procbed:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		bed=rules.make_annotation_bed.output
	output:
		config['output'] + '/odgi/procbed/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_procbed.benchmark.txt'
	shell:
		'''
		odgi procbed  \
		-i {input.og} \
		-b {input.bed} > {output}
		'''

rule annot_names:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.odgi_procbed.output
	output:
		config['output'] + '/odgi/procbed/{region}.names.txt'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{region}.annot_names.benchmark.txt'
	shell:
		'''
		cut -f 4 {input} > {output}
		'''

rule annot_path:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.odgi_procbed.output
	output:
		config['output'] + '/odgi/procbed/{region}.path.txt'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{region}.annot_path.benchmark.txt'
	shell:
		'''
		cut -f 1 {input} | uniq > {output}
		'''

rule odgi_inject:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		bed=rules.odgi_procbed.output
	output:
		config['output'] + '/odgi/inject/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_inject.benchmark.txt'
	shell:
		'''
		odgi inject \
			-i {input.og} \
			-b {input.bed} \
			-o {output}
		'''

rule odgi_flip:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.odgi_inject.output,
		txt=rules.annot_path.output
	output:
		config['output'] + '/odgi/flip/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_flip.benchmark.txt'
	params:

	shell:
		'''
		odgi flip \
			-i {input.og} \
			-o {output} \
			--ref-flips {input.txt}
		'''	

rule odgi_untangle:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.odgi_flip.output,
		txt=rules.annot_names.output
	output:
		config['output'] + '/odgi/untangle/{region}.gggenes.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1736526388'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{region}.odgi_untangle.benchmark.txt'
	shell:
		'''
		odgi untangle \
			-R {input.txt} \
			-i {input.og} \
			-j 0.3 \
			-g > {output}
		'''

rule plot_gggenes:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		tsv=rules.odgi_untangle.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/odgi/untangle/{region}.gggenes.pdf'
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
		'benchmarks/{region}.plot_gggenes.benchmark.txt'
	shell:
		'''
		Rscript \
			workflow/scripts/plotgggenes.r \
			{input.tsv} \
			{input.json} \
			{output}
		'''

rule plot_gggenes_dbscan:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		tsv=rules.odgi_untangle.output,
		json=rules.make_dbscan_clusters.output
	output:
		config['output'] + '/odgi/untangle_dbscan/{region}.gggenes.pdf'
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
		'benchmarks/{region}.plot_gggenes.benchmark.txt'
	shell:
		'''
		Rscript \
			workflow/scripts/plotgggenes.r \
			{input.tsv} \
			{input.json} \
			{output}
		'''