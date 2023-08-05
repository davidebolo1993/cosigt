rule pggb:
    '''
    pggb
    '''
	input:
		'...',
	output:
		'results/pggb/{region}/subgraph.gfa.gz'
	threads:
		config['pggb']['threads']
	container:
		'docker://pangenome/pggb:latest'
	resources:
		mem_mb=config['pggb']['mem_mb'],
		time=config['pggb']['time']
	shell:
		'''
        pggb \
            -i {input} \
            -o {output} \
            -p 90 \
            -s 5k \
            -k 419 \
            -t {threads} \
            -n 200
		'''
