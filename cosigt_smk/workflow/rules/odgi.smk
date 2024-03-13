rule odgi_chop:
	'''
	odgi chop
	'''
	input:
		rules.pggb.output
	output:
		config['output'] + '/odgi/chop/{region}.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	shell:
		'''
		odgi chop  \
		-i {input} \
		-c 32 \
		-o {output}
		'''

rule odgi_view:
	'''
	odgi view
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
	shell:
		'''
		odgi view \
		-i {input} \
		-g > {output}
		'''

rule odgi_paths_matrix:
	'''
	odgi build non-binary haplotype matrix
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
	shell:
		'''
		odgi paths \
		-i {input} \
		-H | cut -f 1,4- | gzip > {output}
		'''

rule odgi_similarity:
	'''
	odgi similarity
	'''
	input:
		rules.pggb.output
	output:
		config['output'] + '/odgi/similarity/{region}.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	shell:
		'''
		odgi similarity \
		-i {input} > {output}
		'''	