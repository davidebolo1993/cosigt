rule pggb:
	'''
	pggb
	'''
	input:
		fasta=rules.odgi_paths_fasta.output,
		index=rules.faidx.output
	output:
		config['output'] + '/pggb/{region}.og'
	threads:
		config['pggb']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	container:
		'docker://pangenome/pggb:latest'
	params:
		prefix=config['output'] + '/pggb/{region}',
		flags=config['pggb']['params']
	shell:
		'''
		pggb \
			-i {input.fasta} \
			-o {params.prefix} \
			-t {threads} \
			{params.flags} \
		&& mv {params.prefix}/*smooth.final.og {output}
		'''