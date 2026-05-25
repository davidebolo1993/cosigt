if LONG_READ_PRESET is not None:

	rule minimap2_reads_index:
		'''
		https://github.com/lh3/minimap2
		- Build a preset-specific minimap2 index for allele contigs
		- Reused across all samples for this region and long-read preset
		'''
		input:
			ref_fasta=rules.bedtools_getfasta.output.fasta
		output:
			outpath(f"minimap2/{READ_MODE_LABEL}/index/{{chr}}/{{region}}/{{region}}.mmi")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['reads']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['minimap2']['reads']['runtime']
		container:
			'docker://davidebolo1993/minimap2:2.31'
		conda:
			'../envs/minimap2.yaml'
		benchmark:
			f'benchmarks/{{chr}}.{{region}}.{READ_MODE_LABEL}.minimap2_reads_index.benchmark.txt'
		params:
			preset=lambda wildcards: long_read_preset()
		shell:
			'''
			minimap2 \
				-x {params.preset} \
				-d {output} \
				{input.ref_fasta}
			'''

	rule minimap2_reads_samtools_sort:
		'''
		https://github.com/lh3/minimap2
		https://github.com/samtools/samtools
		- Re-align long reads to allele contigs with a read-technology preset
		- Sort and convert to indexed CRAM for graph injection
		'''
		input:
			ref_fasta=rules.bedtools_getfasta.output.fasta,
			ref_index=rules.minimap2_reads_index.output,
			sample_fasta=rules.combine_mapped_unmapped.output
		output:
			cram=temp(outpath(f"minimap2/{READ_MODE_LABEL}/{{sample}}/{{chr}}/{{region}}/{{region}}.realigned.cram")),
			crai=temp(outpath(f"minimap2/{READ_MODE_LABEL}/{{sample}}/{{chr}}/{{region}}/{{region}}.realigned.cram.crai"))
		threads:
			config['minimap2']['reads']['threads']
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['reads']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['minimap2']['reads']['runtime']
		container:
			'docker://davidebolo1993/minimap2:2.31'
		conda:
			'../envs/minimap2.yaml'
		benchmark:
			f'benchmarks/{{sample}}.{{chr}}.{{region}}.{READ_MODE_LABEL}.minimap2_reads_samtools_sort.benchmark.txt'
		params:
			tmp_prefix=outpath(f"minimap2/{READ_MODE_LABEL}/{{sample}}/{{chr}}/{{region}}"),
			preset=lambda wildcards: long_read_preset()
		shell:
			'''
			minimap2 \
				-a \
				-x {params.preset} \
				-t {threads} \
				{input.ref_index} \
				{input.sample_fasta} | \
			samtools sort \
				-T {params.tmp_prefix} | \
			samtools view \
				-T {input.ref_fasta} \
				-C \
				-o {output.cram} \
				--write-index
			'''
