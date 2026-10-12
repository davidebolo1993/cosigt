if READ_MODE == 'short':

	rule bwamem3_index:
		'''
		https://github.com/fg-labs/bwa-mem3
		- Build index for the contigs
		- Written under an aligner-specific prefix rather than next to the shared
		  allele FASTA, so that switching read_mode cannot leave behind stale
		  .pac/.ann/.amb files in the incompatible bwa format
		- No .0123: bwa-mem3 fetches bases from .pac and never reads the unpacked
		  reference, so it is only written with --emit-unpacked-ref
		- -t is explicit because the default is every core the node reports
		'''
		input:
			rules.bedtools_getfasta.output.fasta
		output:
			multiext(outpath("bwa-mem3/index/{chr}/{region}/{region}.fasta.gz"), '.bwt.2bit.64', '.pac', '.ann', '.amb')
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem3']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa-mem3']['runtime']
		container:
			'docker://davidebolo1993/bwa-mem3:0.14.1'
		conda:
			'../envs/bwa-mem3.yaml'
		benchmark:
			'benchmarks/{chr}.{region}.bwamem3_index.benchmark.txt'
		params:
			prefix=outpath("bwa-mem3/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			bwa-mem3 index -t {threads} -p {params.prefix} {input}
			'''

	rule bwamem3_mem_samtools_sort:
		'''
		https://github.com/fg-labs/bwa-mem3
		https://github.com/samtools/samtools
		- Re-align original short-reads to the contigs, keeping up to 10k multi-mappings
		- The mapped and the k-mer-filtered unmapped reads are streamed into
		  bwa-mem3 together rather than first concatenated by a job of their own,
		  which was one extra job per (sample, region) for a `cat`
		- Sort alignment
		- Convert to .cram and index at the same time
		- bwa-mem3 is a fork of bwa-mem2 run here in its default mode, which keeps
		  bwa-mem2's alignments apart from upstream correctness fixes, resolves
		  exact ties between equally scoring haplotypes differently, and adds MQ
		  and HN tags. None of that reaches the genotype: gfainject takes every
		  alternative hit from XA and gafpack weights each read across them, so
		  which tied haplotype is primary is irrelevant to node coverage.
		'''
		input:
			ref_fasta=rules.bedtools_getfasta.output.fasta,
			ref_fai=rules.bwamem3_index.output,
			fasta_mapped=rules.samtools_fasta_mapped.output,
			fasta_unmapped=rules.kfilt_filter_unmapped.output
		output:
			cram=temp(outpath("bwa-mem3/{sample}/{chr}/{region}/{region}.realigned.cram")),
			crai=temp(outpath("bwa-mem3/{sample}/{chr}/{region}/{region}.realigned.cram.crai"))
		group:
			"genotype"
		threads:
			config['bwa-mem3']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem3']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['bwa-mem3']['runtime'],
			# bwa-mem2 aborts with "assert failed for seqPair size" when too many
			# seeds share a position. That is the normal case here: the reference
			# is an allele panel in which every locus repeats once per haplotype,
			# so seed occurrence counts scale with panel size. Lowering -c is the
			# documented workaround (bwa-mem2 issue #269). bwa-mem3 is a fork and
			# documents no fix for it, so keep halving -c on each retry. Carried
			# as a resource because only resource functions are passed `attempt`.
			max_occ=lambda wildcards, attempt: max(
				config['bwa-mem3']['min_max_occ'],
				config['bwa-mem3']['max_occ'] // (2 ** (attempt - 1))
			)
		container:
			'docker://davidebolo1993/bwa-mem3:0.14.1'
		conda:
			'../envs/bwa-mem3.yaml'
		benchmark:
			'benchmarks/{sample}.{chr}.{region}.bwamem3_mem_samtools_sort.benchmark.txt'
		retries:
			config['bwa-mem3']['retries']
		params:
			tmp_prefix=outpath("bwa-mem3/{sample}/{chr}/{region}"),
			index_prefix=outpath("bwa-mem3/index/{chr}/{region}/{region}.fasta.gz")
		shell:
			'''
			echo "bwa-mem3 mem -c {resources.max_occ} (seed-occurrence cap, halved on each retry)" >&2
			cat {input.fasta_mapped} {input.fasta_unmapped} | bwa-mem3 mem \
			-t {threads} \
			-p \
			-h 10000 \
			-c {resources.max_occ} \
			{params.index_prefix} \
			/dev/stdin | samtools sort \
			-@ {threads} \
			-T {params.tmp_prefix} | \
			samtools view \
			-T {input.ref_fasta} \
			-C \
			-o {output.cram} \
			--write-index
			'''
