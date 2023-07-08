rule odgi_chop:
	'''
	odgi chop
	'''
	input:
		config['graph']
	output:
		'resources/odgi/z.gfa'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi chop  \
		-i {input} \
		-c 32 \
		-o - | odgi view \
		-i - \
		-g > {output}
		'''

rule odgi_build:
	'''
	odgi build non-binary haplotype matrix
	'''
	input:
		rules.odgi_chop.output
	output:
		'resources/odgi/z.paths.tsv.gz'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi build \
		-g {input} \
		-o - | odgi paths \
		-i - \
		-H | cut -f 1,4- | pigz > {output}
		'''

rule odgi_paths:
	'''
	odgi paths
	'''
	input:
		config['graph']
	output:
		'resources/odgi/z.fa'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi paths \
		-i {input} \
		-f > {output}
		'''

rule odgi_similarity:
	'''
	odgi similarity
	'''
	input:
		config['graph']
	output:
		'resources/odgi/hapdiff.tsv'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi build \
		-g {input} \
		-o - | odgi similarity \
		-i - > {output}
		'''	