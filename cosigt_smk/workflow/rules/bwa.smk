rule bwa_index:
	'''
	https://github.com/lh3/bwa
	- Build index for the contigs
	'''
	input:
		rules.bedtools_getfasta.output.fasta
	output:
		multiext(config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz', '.bwt', '.pac', '.ann', '.amb', '.sa')	
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa']['time']
	container:
		'docker://davidebolo1993/bwa:0.7.18'
	conda:
		'../envs/bwa.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.bwa_index.benchmark.txt'
	shell:
		'''
		bwa index {input}
		'''

rule bwa_aln:
	'''
	https://github.com/lh3/bwa
	- Re-align original short-reads to the contigs, using parameters for ancient DNA
	'''
	input:
		ref_fasta=rules.bedtools_getfasta.output.fasta,
		ref_fai=rules.bwa_index.output,
		sample_fasta=rules.samtools_fasta_mapped.output
	output:
		config['output'] + '/bwa/{sample}/{chr}/{region}/{region}.realigned.sai'
	threads:
		config['bwa']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa']['time']
	container:
		'docker://davidebolo1993/bwa:0.7.18'
	conda:
		'../envs/bwa.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.bwa_aln_samse_samtools_sort.benchmark.txt'
	shell:
		'''
		bwa aln \
		-t {threads} \
		-0 \
		-l 1024 \
		-n 0.01 \
		-o 2 \
		{input.ref_fasta} \
		{input.sample_fasta} > {output}
		'''

rule bwa_samse_samtools_sort:
	'''
	https://github.com/lh3/bwa
	https://github.com/samtools/samtools
	- Sai to bam
	- Bam to cram and indexing
	'''
	input:
		sai=rules.bwa_aln.output,
		ref_fasta=rules.bedtools_getfasta.output.fasta,
		ref_fai=rules.bwa_index.output,
		sample_fasta=rules.samtools_fasta_mapped.output
	output:
		config['output'] + '/bwa/{sample}/{chr}/{region}/{region}.realigned.cram'
	threads:
		config['samtools']['fasta_mapped']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['time']
	container:
		'docker://davidebolo1993/bwa:0.7.18'
	conda:
		'../envs/bwa.yaml'	
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.bwa_samse_samtools_sort.benchmark.txt'
	params:
		tmp_prefix=config['output'] + '/bwa/{sample}/{chr}/{region}'
	shell:
		'''
		bwa samse \
		-n 10000 \
		{input.ref_fasta} \
		{input.sai} \
		{input.sample_fasta} | \
		samtools sort \
		-T {params.tmp_prefix} | \
		samtools view \
		-T {input.ref_fasta} \
		-C \
		-o {output} \
		--write-index
		'''

rule wally_viz_alignments:
	'''
	https://github.com/tobiasrausch/wally
	- Visualise reads-to-contig alignment for each contig
	'''
	input:
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai,
		sample=rules.bwa_samse_samtools_sort.output
	output:
		config['output'] + '/wally/{sample}/{chr}/{region}/wally.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['small']['time']
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
		cut -f 1,2 {input.fai} | awk '{{print $1,1,$2}}' OFS="\\t" > plot.bed
		wally \
			region \
			-g {input.fasta} \
			-R plot.bed \
			-q 0 \
			-c \
			-u \
			-p \
			-y 512 \
			-x 4096 \
			{input.sample}
		rm plot.bed
		touch {output}
		'''	
