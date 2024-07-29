rule megadepth_bam_to_bigwig:
	'''
	https://github.com/ChristopherWilks/megadepth
	'''
	input:
		rules.bwamem2_mem_samtools_sort.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{region}.realigned.bam.all.bw'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{sample}.{region}.megadepth_bam_to_bigwig.benchmark.txt'
	shell:
		'''
		megadepth \
		--bigwig \
		{input}
		'''