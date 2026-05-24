if READ_MODE == 'short':

	rule bwamem2_index:
		'''
		https://github.com/bwa-mem2/bwa-mem2
		- Build index for the contigs
		'''
		input:
			rules.bedtools_getfasta.output.fasta
		output:
			multiext(outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz"), '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa-mem2']['runtime']
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
		- Re-align original short-reads to the contigs, keeping up to 10k multi-mappings
		- Sort alignment
		- Convert to .cram and index at the same time
		'''
		input:
			ref_fasta=rules.bedtools_getfasta.output.fasta,
			ref_fai=rules.bwamem2_index.output,
			sample_fasta=rules.combine_mapped_unmapped.output
		output:
			cram=temp(outpath("bwa-mem2/{sample}/{chr}/{region}/{region}.realigned.cram")),
			crai=temp(outpath("bwa-mem2/{sample}/{chr}/{region}/{region}.realigned.cram.crai"))
		threads:
			config['bwa-mem2']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa-mem2']['runtime']
		container:
			'docker://davidebolo1993/bwa-mem2:2.2.1'
		conda:
			'../envs/bwa-mem2.yaml'
		benchmark:
			'benchmarks/{sample}.{chr}.{region}.bwamem2_mem_samtools_sort.benchmark.txt'
		params:
			tmp_prefix=outpath("bwa-mem2/{sample}/{chr}/{region}")
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
			-o {output.cram} \
			--write-index
			'''
