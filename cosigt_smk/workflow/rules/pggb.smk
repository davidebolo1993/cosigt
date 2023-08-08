rule pggb:
    '''
    pggb
    '''
	input:
		fasta=rules.pgrtk_get_seq.output,
		index=rules.faidx.output
	output:
		'results/pggb/{region}.og'
	threads:
		config['pggb']['threads']
	container:
		'docker://pangenome/pggb:latest'
	resources:
		mem_mb=config['pggb']['mem_mb'],
		time=config['pggb']['time']
	params:
		prefix='results/pggb/{region}'
	shell:
		'''
        pggb \
            -i {input.fasta} \
            -o {params.prefix} \
            -p 90 \
            -s 5k \
            -k 419 \
            -t {threads} \
            -n 200 \
		&& mv {params.prefix}/*smooth.final.og {output}
		'''
