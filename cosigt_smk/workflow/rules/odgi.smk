rule odgi_chop:
	'''
	https://github.com/pangenome/odgi
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
		'docker://pangenome/odgi:1707818641'
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
		rules.odgi_chop.output
	output:
		config['output'] + '/odgi/view/{region}.gfa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
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
		rules.odgi_chop.output
	output:
		config['output'] + '/odgi/paths/matrix/{region}.tsv.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	benchmark:
		'benchmarks/{region}.odgi_paths_matrix.benchmark.txt'
	shell:
		'''
		odgi paths \
		-i {input} \
		-H | \
		cut -f 1,4- | gzip > {output}
		'''

rule odgi_similarity:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/odgi/similarity/{region}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	benchmark:
		'benchmarks/{region}.odgi_similarity.benchmark.txt'
	shell:
		'''
		odgi similarity \
		-i {input} > {output}
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
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{region}.make_clusters.benchmark.txt'
	params:
		prefix=config['output'] + '/cluster/{region}'
	shell:
		'''
		python workflow/scripts/cluster.py \
			{input} \
			{params.prefix}
		'''