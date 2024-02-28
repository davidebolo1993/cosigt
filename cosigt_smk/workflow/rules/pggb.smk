rule pggb:
	'''
	pggb
	'''
	input:
		fasta=rules.odgi_paths_fasta.output,
		index=rules.faidx.output
	output:
		'results/pggb/{region}.og'
	threads:
		config['pggb']['threads']
	container:
		'docker://pangenome/pggb:latest'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	params:
		prefix='results/pggb/{region}'
	shell:
		'''
		pggb \
			-i {input.fasta} \
			-o {params.prefix} \
			-t {threads} \
			-n 94 \
			-c 2 \
		&& mv {params.prefix}/*smooth.final.og {output}
		'''
