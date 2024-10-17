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
	conda:
		'../envs/megadepth.yaml'
	benchmark:
		'benchmarks/{sample}.{region}.megadepth_bam_to_bigwig.benchmark.txt'
	shell:
		'''
		megadepth \
		--bigwig \
		{input}
		'''

rule plot_bigwig_coverage:
	'''
	https://github.com/davidebolo1993/cosigt
	'''	
	input:
		rules.megadepth_bam_to_bigwig.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{region}.realigned.bam.all.pdf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	conda:
		'../envs/plot.yaml'
	shell:
		'''
		Rscript \
		workflow/scripts/plotcoverage.r \
		{input} \
		{output}
		'''