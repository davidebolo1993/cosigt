rule pggb_construct:
	'''
	https://github.com/pangenome/pggb
	'''
	input:
		fasta=rules.samtools_faidx_extract.output,
		fai=rules.samtools_faidx_index.output
	output:
		config['output'] + '/pggb/{region}.og'
	threads:
		config['pggb']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	container:
		'docker://ghcr.io/pangenome/pggb:202410231852260e9c9e'
	conda:
		'../envs/pggb.yaml'
	benchmark:
		'benchmarks/{region}.pggb_construct.benchmark.txt'
	params:
		prefix=config['output'] + '/pggb/{region}',
		flags=config['pggb']['params'],
		tmpdir = config['pggb']['tmpdir']
	shell:
		'''
		pggb \
			-i {input.fasta} \
			-o {params.prefix} \
			-t {threads} \
			-D {params.tmpdir} \
			-n $(wc -l {input.fai}) \
			{params.flags} \
		&& mv {params.prefix}/*smooth.final.og {output}
		'''
