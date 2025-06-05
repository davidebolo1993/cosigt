rule pggb_construct:
	'''
	https://github.com/pangenome/pggb
	'''
	input:
		fasta=rules.samtools_faidx_extract.output,
		fai=rules.samtools_faidx_index.output
	output:
		config['output'] + '/pggb/{chr}/{region}.og'
	threads:
		config['pggb']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['pggb']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['pggb']['time']
	container:
		'docker://ghcr.io/pangenome/pggb:20250423145743e25486'
	conda:
		'../envs/pggb.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pggb_construct.benchmark.txt'
	params:
		prefix=config['output'] + '/pggb/{chr}/{region}',
		flags=config['pggb']['params'],
		tmpdir = config['pggb']['tmpdir'] + '/{chr}/{region}'
	shell:
		'''
		mkdir -p {params.tmpdir}
		pggb \
			-i {input.fasta} \
			-o {params.prefix} \
			-t {threads} \
			-D {params.tmpdir} \
			{params.flags} \
		&& mv {params.prefix}/*smooth.final.og {output}
		rm -rf {params.tmpdir}
		'''
