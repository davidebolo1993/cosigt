rule wally_viz_alignments:
	'''
	https://github.com/tobiasrausch/wally
	- Visualise reads-to-contig alignment for each contig
	'''
	input:
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai,
		sample=realigned_alignment_path
	output:
		outpath("wally/{sample}/{chr}/{region}/wally.done")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://trausch/wally:v0.7.1'
	conda:
		'../envs/wally.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.wally_viz_alignments.benchmark.txt'
	params:
		outdir=outpath("wally/{sample}/{chr}/{region}")
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
