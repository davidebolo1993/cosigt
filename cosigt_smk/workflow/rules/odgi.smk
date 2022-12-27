rule odgi_chop:
	'''
	Odgi chop
	'''
	input:
		config['graph']
	output:
		"resources/odgi/z.gfa"
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
	Odgi build non-binary haplotype matrix
	'''
	input:
		rules.odgi_chop.output
	output:
		"resources/odgi/z.paths.tsv.gz"
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
	Odgi paths
	'''
	input:
		config['graph']
	output:
		"resources/odgi/z.fa"
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
