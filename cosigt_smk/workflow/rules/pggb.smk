rule pggb_construct:
	'''
	https://github.com/pangenome/pggb
	'''
	input:
		fasta=rules.pyfaidx_extract.output,
		index=rules.samtools_faidx_index.output
	output:
		config['output'] + '/pggb/{region}.og'
	threads:
		config['pggb']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	container:
		'docker://pangenome/pggb:latest'
	benchmark:
		'benchmarks/{region}.pggb_construct.benchmark.txt'
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