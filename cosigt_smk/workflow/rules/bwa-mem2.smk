if READ_MODE == 'short':

	rule bwamem2_index:
		'''
		https://github.com/bwa-mem2/bwa-mem2
		- Build index for the contigs
		- Written under an aligner-specific prefix rather than next to the shared
		  allele FASTA, so that switching read_mode cannot leave behind stale
		  .pac/.ann/.amb files in the incompatible bwa format
		'''
		input:
			rules.bedtools_getfasta.output.fasta
		output:
			multiext(outpath("bwa-mem2/index/{chr}/{region}/{region}.fasta.gz"), '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')
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
		params:
			prefix=outpath("bwa-mem2/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			bwa-mem2 index -p {params.prefix} {input}
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
		group:
			"genotype"
		threads:
			config['bwa-mem2']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa-mem2']['runtime'],
			# bwa-mem2 aborts with "assert failed for seqPair size" when too many
			# seeds share a position. That is the normal case here: the reference
			# is an allele panel in which every locus repeats once per haplotype,
			# so seed occurrence counts scale with panel size. Lowering -c is the
			# documented workaround (bwa-mem2 issue #269), so halve it on each
			# retry rather than editing this rule by hand. Carried as a resource
			# because only resource functions are passed `attempt`.
			max_occ=lambda wildcards, attempt: max(
				config['bwa-mem2']['min_max_occ'],
				config['bwa-mem2']['max_occ'] // (2 ** (attempt - 1))
			)
		container:
			'docker://davidebolo1993/bwa-mem2:2.2.1'
		conda:
			'../envs/bwa-mem2.yaml'
		benchmark:
			'benchmarks/{sample}.{chr}.{region}.bwamem2_mem_samtools_sort.benchmark.txt'
		retries:
			config['bwa-mem2']['retries']
		params:
			tmp_prefix=outpath("bwa-mem2/{sample}/{chr}/{region}"),
			index_prefix=outpath("bwa-mem2/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			echo "bwa-mem2 mem -c {resources.max_occ} (seed-occurrence cap, halved on each retry)" >&2
			bwa-mem2 mem \
			-t {threads} \
			-p \
			-h 10000 \
			-c {resources.max_occ} \
			{params.index_prefix} \
			{input.sample_fasta} | samtools sort \
			-@ {threads} \
			-T {params.tmp_prefix} | \
			samtools view \
			-T {input.ref_fasta} \
			-C \
			-o {output.cram} \
			--write-index
			'''
