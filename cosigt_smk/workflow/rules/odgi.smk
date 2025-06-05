rule odgi_view:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/view/{chr}/{region}.gfa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_view.benchmark.txt'
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
		config['output'] + '/odgi/paths/matrix/{chr}/{region}.tsv.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_paths_matrix.benchmark.txt'
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
		config['output'] + '/odgi/view/{chr}/{region}.node.length.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.get_node_length.benchmark.txt'
	shell:
		'''
		grep '^S' {input} | \
		awk '{{print("node."$2,length($3))}}' OFS="\\t" > {output}
		'''

checkpoint generate_submasks:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.get_nodes_length.output
	output:
		directory(config['output'] + '/odgi/paths/matrix/{chr}/{region}_submasks')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	benchmark:
		'benchmarks/{chr}.{region}.generate_submasks.benchmark.txt'
	shell:
		'''
		bash workflow/scripts/make_mask_files.sh {input} {output} 10000
		'''

def get_submasks(wildcards):
	'''
	extract sub-masks
	'''
	region=wildcards.region
	checkpoint_output = checkpoints.generate_submasks.get(region=region).output[0]
	submask_files = glob(checkpoint_output + '/*tsv')
	return [os.path.basename(f).split('.')[0] for f in submask_files]


rule odgi_dissimilarity_submasks:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		mask=config['output'] + '/odgi/paths/matrix/{chr}/{region}_submasks/{num}.tsv'
	output:
		config['output'] + '/odgi/dissimilarity/{chr}/{region}_submasks/{num}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.{num}.odgi_dissimilarity_submasks.benchmark.txt'
	shell:
		'''
		odgi similarity \
		-i {input.og} \
		-m {input.mask} \
		--all \
		--distances > {output}
		'''

rule collect_submasks_dissimilarities:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/odgi/dissimilarity{chr}/{region}_submasks/{num}.tsv', chr=wildcards.chr, region=wildcards.region, num=get_submasks(wildcards))
	output:
		config['output'] + '/odgi/dissimilarity/{chr}/{region}/submasks.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.collect_submasks_dissimilarities.benchmark.txt'
	shell:
		'''
		touch {output}
		'''

rule filter_nodes:
	'''
	https://github.com/davidebolo1993/cosigt 
	'''	
	input:
		graph_cov=rules.odgi_paths_matrix.output,
		node_length=rules.get_nodes_length.output
	output:
		config['output'] + '/odgi/paths/matrix/{chr}/{region}.mask.tsv'
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
		{output}
		'''	

rule odgi_dissimilarity:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		mask=rules.filter_nodes.output
	output:
		config['output'] + '/odgi/dissimilarity/{chr}/{region}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
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
		--distances > {output}
		'''

rule make_clusters:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		matrix=rules.odgi_dissimilarity.output,
		filter=rules.filter_nodes.output
	output:
		config['output'] + '/cluster/{chr}/{region}.clusters.json'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_clusters.benchmark.txt'
	params:
		threshold_file=config['output'] + '/odgi/paths/matrix/{chr}/{region}.shared.tsv'
	shell:
		'''
		Rscript workflow/scripts/cluster.r \
			{input.matrix} \
			{output} \
			automatic \
			$(cut -f 3 {params.threshold_file} | tail -1)
		'''

rule viz_odgi:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		graph_cov=rules.graph_cov.output,
		json=rules.make_clusters.output,
		nodes_length=rules.get_nodes_length.output
	output:
		config['output'] + '/odgi/viz/{chr}/{region}.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.viz_odgi.benchmark.txt'
	shell:
		'''
		Rscript \
			workflow/scripts/viz_odgi.r \
			{input.graph_cov} \
			{input.json} \
			{input.nodes_length} \
			{output}
		'''

rule subset_gtf:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		gtf=config['annotations'],
		bed=rules.make_reference_bed.output
	output:
		config['output'] + '/annotations/{chr}/{region}.gtf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.subset_gtf.benchmark.txt'
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
		config['output'] + '/annotations/{chr}/{region}.bed'
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
		'benchmarks/{chr}.{region}.make_annotation_bed.benchmark.txt'
	params:
		pansn=config['pansn_prefix'] + '{chr}'
	shell:
		'''
		Rscript workflow/scripts/annotate.r \
		{input} \
		{params.pansn} \
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
		config['output'] + '/odgi/procbed/{chr}/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_procbed.benchmark.txt'
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
		config['output'] + '/odgi/procbed/{chr}/{region}.names.txt'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.annot_names.benchmark.txt'
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
		config['output'] + '/odgi/procbed/{chr}/{region}.path.txt'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.annot_path.benchmark.txt'
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
		config['output'] + '/odgi/inject/{chr}/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_inject.benchmark.txt'
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
		config['output'] + '/odgi/flip/{chr}/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_flip.benchmark.txt'
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
		config['output'] + '/odgi/untangle/{chr}/{region}.gggenes.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_untangle.benchmark.txt'
	shell:
		'''
		odgi untangle \
			-R {input.txt} \
			-i {input.og} \
			-j 0.3 \
			-g \
			-m 1000 \
			-e 1000 > {output}
		'''

rule plot_gggenes:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		tsv=rules.odgi_untangle.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/odgi/untangle/{chr}/{region}.gggenes.pdf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.plot_gggenes.benchmark.txt'
	shell:
		'''
		Rscript \
			workflow/scripts/plotgggenes.r \
			{input.tsv} \
			{input.json} \
			{output}
		'''
