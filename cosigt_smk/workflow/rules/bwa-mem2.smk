rule bwamem2_index:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	'''
	input:
		rules.samtools_faidx_extract.output
	output:
		multiext(config['output'] + '/samtools/faidx/{chr}/{region}.fasta', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')	
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/bwa-mem2:2.2.1'
	conda:
		'../envs/bwa-mem2.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.bwamem2_index.benchmark.txt'
	shell:
		'''
		bwa-mem2 index {input}
		'''

rule bwamem2_mem_samtools_sort:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	https://github.com/samtools/samtools
	'''
	input:
		ref_fasta=rules.samtools_faidx_extract.output,
		ref_fai=rules.bwamem2_index.output,
		sample_fasta=rules.samtools_fasta_mapped.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{chr}/{region}.realigned.bam'
	threads:
		config['bwa-mem2']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/bwa-mem2:2.2.1'
	conda:
		'../envs/bwa-mem2.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.bwamem2_mem_samtools_sort.benchmark.txt'
	params:
		tmp_prefix=config['output'] + '/bwa-mem2/{sample}/{chr}/{region}'
	shell:
		'''
		bwa-mem2 mem \
		-t {threads} \
		-p \
		-h 10000 \
		{input.ref_fasta} \
		{input.sample_fasta} | samtools sort \
		-T {params.tmp_prefix} > {output}
		'''

rule samtools_index_realigned_bam:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.bwamem2_mem_samtools_sort.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{chr}/{region}.realigned.bam.bai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.samtools_index_realigned_bam.benchmark.txt'
	shell:
		'''
		samtools index {input}
		'''	

rule make_wally_plot_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.samtools_faidx_index.output
	output:
		config['output'] + '/samtools/faidx/{chr}/{region}.plot.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.make_wally_plot_bed.benchmark.txt'
	shell:
		'''
		cut -f 1,2 {input} | awk '{{print $1,1,$2}}' OFS="\\t" > {output}
		'''
	
rule wally_viz_alignments:
	'''
	https://github.com/tobiasrausch/wally
	'''
	input:
		bed=rules.make_wally_plot_bed.output,
		fasta=rules.samtools_faidx_extract.output,
		sample=rules.bwamem2_mem_samtools_sort.output,
		sample_index=rules.samtools_index_realigned_bam.output
	output:
		config['output'] + '/wally/{sample}/{chr}/{region}/wally.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://trausch/wally:v0.7.1'
	conda:
		'../envs/wally.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.wally_viz_alignments.benchmark.txt'
	params:
		outdir=config['output'] + '/wally/{sample}/{chr}/{region}'
	shell:
		'''
		cd {params.outdir}
		wally \
			region \
			-g {input.fasta} \
			-R {input.bed} \
			-q 0 \
			-c \
			-u \
			-p \
			-y 512 \
			-x 4096 \
			{input.sample}
		touch {output}
		'''	
