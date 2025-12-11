rule bwamem2_index:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	- Build index for the contigs
	'''
	input:
		rules.copy_fasta_over.output.fasta
	output:
		temp(multiext(config['output'] + '/alleles/{chr}/{region}/{region}.fasta.gz', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123'))
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
	params:
		outdir=config['output']+ '/alleles/{chr}/{region}'
	shell:
		'''
		bwa-mem2 index {input}
		'''
	
rule bwamem2_mem_samtools_sort:
	'''
	https://github.com/bwa-mem2/bwa-mem2
	https://github.com/samtools/samtools
	- Re-align original short-reads to the contigs, keeping up to 10k multi-mappings
	- Sort alignment
	- Convert to .cram and index at the same time
	'''
	input:
		ref_fasta=rules.copy_fasta_over.output.fasta,
		ref_fai=rules.bwamem2_index.output,
		sample_fasta=rules.combine_mapped_unmapped.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{chr}/{region}/{region}.realigned.cram'
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
		fasta=rules.copy_fasta_over.output.fasta,
		fai=rules.copy_fasta_over.output.fai,
		sample=rules.bwamem2_mem_samtools_sort.output
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
