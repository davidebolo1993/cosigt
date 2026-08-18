if READ_MODE == 'ancient':

	rule bwa_index:
		'''
		https://github.com/lh3/bwa
		- Build index for allele contigs for ancient-DNA realignment
		- Written under an aligner-specific prefix rather than next to the shared
		  allele FASTA, so that switching read_mode cannot leave behind stale
		  .pac/.ann/.amb files in the incompatible bwa-mem2 format
		'''
		input:
			rules.bedtools_getfasta.output.fasta
		output:
			multiext(outpath("bwa/index/{chr}/{region}/{region}.fasta.gz"), '.bwt', '.pac', '.ann', '.amb', '.sa')
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa']['runtime']
		container:
			'docker://davidebolo1993/bwa:0.7.18'
		conda:
			'../envs/bwa.yaml'
		benchmark:
			'benchmarks/{chr}.{region}.bwa_index.benchmark.txt'
		params:
			prefix=outpath("bwa/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			bwa index -p {params.prefix} {input}
			'''


	rule bwa_aln:
		'''
		https://github.com/lh3/bwa
		- Re-align original short reads to allele contigs using ancient-DNA parameters
		'''
		input:
			ref_fasta=rules.bedtools_getfasta.output.fasta,
			ref_fai=rules.bwa_index.output,
			sample_fasta=rules.samtools_fasta_mapped.output
		output:
			temp(outpath("bwa/{sample}/{chr}/{region}/{region}.realigned.sai"))
		group:
			"genotype"
		threads:
			config['bwa']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa']['runtime']
		container:
			'docker://davidebolo1993/bwa:0.7.18'
		conda:
			'../envs/bwa.yaml'
		benchmark:
			'benchmarks/{sample}.{chr}.{region}.bwa_aln.benchmark.txt'
		params:
			index_prefix=outpath("bwa/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			bwa aln \
				-t {threads} \
				-0 \
				-l 1024 \
				-n 0.01 \
				-o 2 \
				{params.index_prefix} \
				{input.sample_fasta} > {output}
			'''


	rule bwa_samse_samtools_sort:
		'''
		https://github.com/lh3/bwa
		https://github.com/samtools/samtools
		- Convert ancient-DNA BWA sai output to sorted/indexed CRAM
		'''
		input:
			sai=rules.bwa_aln.output,
			ref_fasta=rules.bedtools_getfasta.output.fasta,
			ref_fai=rules.bwa_index.output,
			sample_fasta=rules.samtools_fasta_mapped.output
		output:
			cram=temp(outpath("bwa/{sample}/{chr}/{region}/{region}.realigned.cram")),
			crai=temp(outpath("bwa/{sample}/{chr}/{region}/{region}.realigned.cram.crai"))
		group:
			"genotype"
		threads:
			config['samtools']['fasta_mapped']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['runtime']
		container:
			'docker://davidebolo1993/bwa:0.7.18'
		conda:
			'../envs/bwa.yaml'
		benchmark:
			'benchmarks/{sample}.{chr}.{region}.bwa_samse_samtools_sort.benchmark.txt'
		params:
			tmp_prefix=outpath("bwa/{sample}/{chr}/{region}"),
			index_prefix=outpath("bwa/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			bwa samse \
				-n 10000 \
				{params.index_prefix} \
				{input.sai} \
				{input.sample_fasta} | \
			samtools sort \
				-@ {threads} \
				-T {params.tmp_prefix} | \
			samtools view \
				-T {input.ref_fasta} \
				-C \
				-o {output.cram} \
				--write-index
			'''
