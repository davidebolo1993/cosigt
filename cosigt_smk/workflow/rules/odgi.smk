import os

rule odgi_chop:
	'''
	odgi chop
	'''
	input:
		rules.pggb.output
	output:
		'results/odgi/chop/{region}.og'
	threads:
		1
	container:
		'docker://pangenome/odgi:1689324294'
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
		'results/odgi/view/{region}.gfa'
	threads:
		1
	container:
		'docker://pangenome/odgi:1689324294'
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
		'results/odgi/paths/matrix/{region}.tsv.gz'
	threads:
		1
	container:
		'docker://pangenome/odgi:1689324294'
	shell:
		'''
		odgi paths \
		-i {input} \
		-H | cut -f 1,4- | pigz > {output}
		'''

rule odgi_paths:
	'''
	odgi paths
	'''
	input:
		rules.pggb.output
	output:
		'results/odgi/paths/fasta/{region}.fasta'
	threads:
		1
	container:
		'docker://pangenome/odgi:1689324294'
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
		rules.pggb.output
	output:
		'results/odgi/similarity/{region}.tsv'
	threads:
		1
	container:
		'docker://pangenome/odgi:1689324294'
	shell:
		'''
		odgi similarity \
		-i {input} > {output}
		'''	