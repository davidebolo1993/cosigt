#USED ONLY FOR BENCHMARKING - DO NOT USE OTHERWISE
#
#Each sample's own haplotypes are part of the pangenome, so the haplotypes cosigt
#predicts can be compared directly against them. Sequences are compared with
#edlib and reported as QV and error rate; there is no clustering-agreement step.

def benchmark_graph(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- leave_zero_out: the genotyping graph already contains the truth
	- leave_all_out: the truth is absent from the genotyping graph by design and
	  comes from a pre-built per-region graph listed in the truth_graphs table
	'''
	if BENCHMARK_MODE == 'leave_all_out':
		return truth_graph_path(wildcards)
	return outpath("pggb", wildcards.chr, wildcards.region, f"{wildcards.region}.og")


rule odgi_flip_pggb_graph_to_fasta:
	'''
	https://github.com/pangenome/odgi
	- Orient the haplotypes with respect to the target
	- Output fasta file and its index
	- Orientation matters because the edit distance between a sequence and the
	  reverse complement of its match is meaningless
	- Both predicted and true haplotypes are taken from this one FASTA, so they
	  are guaranteed to share an orientation frame
	'''
	input:
		benchmark_graph
	output:
		fasta=outpath("benchmark/{chr}/{region}/{region}.flipped.fasta"),
		fai=outpath("benchmark/{chr}/{region}/{region}.flipped.fasta.fai")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://pangenome/odgi:1753347183'
	conda:
		'../envs/odgi.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.odgi_flip_pggb_graph_to_fasta.benchmark.txt'
	params:
		pansn=config['pansn_prefix'],
		refpath=outpath("benchmark/{chr}/{region}/ref_path.txt")
	shell:
		'''
		odgi paths \
			-i {input} \
			-L | grep {params.pansn} > {params.refpath}
		odgi flip \
			-i {input} \
			-o - \
			--ref-flips {params.refpath} | \
		odgi paths \
			-i - \
			-f | sed 's/_inv$//g' > {output.fasta}
		rm {params.refpath}
		samtools faidx {output.fasta}
		'''


rule benchmark_prepare_qv:
	'''
	https://github.com/davidebolo1993/cosigt
	- Collect, per sample, the predicted and the true haplotype sequences
	- One job per region, looping over samples internally
	'''
	input:
		fasta=rules.odgi_flip_pggb_graph_to_fasta.output.fasta,
		fai=rules.odgi_flip_pggb_graph_to_fasta.output.fai,
		# Names of the haplotypes the genotyping graph actually offered. In
		# leave_all_out these are a strict subset of the FASTA above, and they
		# define the candidate set the oracle is allowed to pick from.
		panel=rules.bedtools_getfasta.output.fai,
		genotypes=lambda wildcards: expand(
			outpath("cosigt/{sample}/{chr}/{region}/{region}.cosigt_genotype.tsv"),
			sample=config['samples'],
			chr=wildcards.chr,
			region=wildcards.region
		)
	output:
		manifest=temp(outpath("benchmark/{chr}/{region}/qv/manifest.tsv")),
		sequences=temp(directory(outpath("benchmark/{chr}/{region}/qv/sequences")))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/samtools:1.23.1'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.benchmark_prepare_qv.benchmark.txt'
	params:
		outdir=outpath("benchmark/{chr}/{region}/qv"),
		mode=BENCHMARK_MODE
	shell:
		'''
		mkdir -p {output.sequences}
		bash workflow/scripts/benchmark_prepare.sh \
			{input.fasta} \
			{input.panel} \
			{params.mode} \
			{params.outdir} \
			{input.genotypes}
		'''


rule benchmark_qv:
	'''
	https://github.com/Martinsos/edlib
	https://github.com/davidebolo1993/cosigt
	- Score predicted against true haplotypes with edlib
	- Report QV and error rate for the better of the two assignments
	'''
	input:
		manifest=rules.benchmark_prepare_qv.output.manifest,
		sequences=rules.benchmark_prepare_qv.output.sequences
	output:
		outpath("benchmark/{chr}/{region}/{region}.qv.tsv")
	threads:
		config['benchmark']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['benchmark']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['benchmark']['runtime']
	container:
		'docker://davidebolo1993/edlib:1.2.7'
	benchmark:
		'benchmarks/{chr}.{region}.benchmark_qv.benchmark.txt'
	params:
		indir=outpath("benchmark/{chr}/{region}/qv"),
		region='{region}',
		gene=lambda wildcards: region_annotation(wildcards.region),
		mode=BENCHMARK_MODE
	shell:
		'''
		bash workflow/scripts/benchmark_qv.sh \
			{params.indir} \
			{params.region} \
			{params.gene:q} \
			{params.mode} \
			{threads} \
			{output}
		'''


def get_all_qv_tables(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- Per-region QV tables for every configured region
	'''
	return [
		outpath("benchmark", REGION_ROWS[region]["chrom"], region, f"{region}.qv.tsv")
		for region in REGION_ORDER
	]


rule benchmark_table:
	'''
	https://github.com/davidebolo1993/cosigt
	- Concatenate the per-region QV tables into a single table
	'''
	input:
		get_all_qv_tables
	output:
		outpath("benchmark/benchmark.qv.tsv")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	benchmark:
		'benchmarks/benchmark_table.benchmark.txt'
	shell:
		'''
		head -n 1 {input[0]} > {output}
		for file in {input}; do
			tail -n +2 "$file" >> {output}
		done
		'''


rule plot_benchmark:
	'''
	https://github.com/davidebolo1993/cosigt
	- Summarise QV per region
	'''
	input:
		rules.benchmark_table.output
	output:
		outpath("benchmark/benchmark.qv.png")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/plot_benchmark.benchmark.txt'
	params:
		mode=BENCHMARK_MODE,
		max_bars=config['benchmark']['max_bars_per_row']
	shell:
		'''
		Rscript \
			workflow/scripts/plot_qv.r \
			{input} \
			{output} \
			{params.mode} \
			{params.max_bars}
		'''
