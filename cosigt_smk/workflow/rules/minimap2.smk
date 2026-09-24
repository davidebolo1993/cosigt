rule pansnspec_target:
	'''
	https://github.com/samtools/samtools
	- Extract the reference chromosome from the reference genome file
	- Convert chromosome name, adapting to PanSN specification (sample#haplotype#contig)
	- Compress with bgzip
	- Build index
	'''
	input:
		config['reference']
	output:
		fasta=outpath("minimap2/{chr}/{chr}.fasta.gz"),
		fai=outpath("minimap2/{chr}/{chr}.fasta.gz.fai")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://davidebolo1993/samtools:1.23.1'
	benchmark:
		'benchmarks/{chr}.pansnspec_target.benchmark.txt'
	conda:
		'../envs/samtools.yaml'
	params:
		pansn=config['pansn_prefix'] + '{chr}'
	shell:
		'''
		samtools faidx \
		{input} \
		{wildcards.chr} | \
		sed "1 s/^.*$/>{params.pansn}/" | \
		bgzip -c > {output.fasta}
		samtools faidx {output.fasta}
		'''

rule minimap2_align_batches:
	'''
	https://github.com/samtools/samtools
	https://github.com/lh3/minimap2
	- Extract the contigs of one PanSN sample (a batch) from the original assemblies
	- Align them (queries) to the target (reference chromosome), streaming the
	  extracted contigs into minimap2 rather than writing them out first
	- Compress with bgzip
	- The batches are read from the assembly .fai at parse time (see
	  assembly_batches), so the contig list is written here rather than by a
	  checkpoint. A checkpoint re-runs the post-processing of the whole DAG each
	  time one completes, which at cohort scale costs as much as building it.
	- Extraction used to be a job of its own. It is a cheap step that is simply
	  redone if the alignment fails, and folding it in removes a job per batch
	  and halves the batch jobs that every group check walks through.
	'''
	input:
		target_fasta=rules.pansnspec_target.output.fasta,
		target_fai=rules.pansnspec_target.output.fai,
		queries_fasta=assembly_fasta_path,
		queries_fai=assembly_fai_path
	output:
		paf=temp(outpath("minimap2/{chr}/batches/paf/{batch}.paf.gz")),
		ids=temp(outpath("minimap2/{chr}/batches/ids/{batch}.txt"))
	threads:
		config['minimap2']['avo']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['avo']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['minimap2']['avo']['runtime']
	container:
		'docker://davidebolo1993/minimap2:2.31'
	conda:
		'../envs/minimap2.yaml'
	benchmark:
		'benchmarks/{chr}.{batch}.minimap2_align_batches.benchmark.txt'
	params:
		batch='{batch}'
	shell:
		'''
		awk -F '\\t' -v batch={params.batch:q} '{{split($1, parts, "#")}} parts[1] == batch {{print $1}}' {input.queries_fai} > {output.ids}
		samtools faidx -r {output.ids} {input.queries_fasta} | \
		minimap2 \
			-x asm20 \
			--eqx \
			-c \
			-t {threads} \
			{input.target_fasta} \
			- | bgzip -c > {output.paf}
		'''

def get_paf_files(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- One PAF per PanSN sample in this chromosome's assemblies
	'''
	return expand(
		outpath("minimap2/{chr}/batches/paf/{batch}.paf.gz"),
		chr=wildcards.chr,
		batch=assembly_batches(wildcards.chr)
	)

rule merge_paf_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	- Concatenate the paf files for each chromosome together
	- Ensure temp files are cleaned up after merging
	- Index .paf - since this is required by impg downstream
	'''
	input:
		get_paf_files
	output:
		paf=outpath("minimap2/{chr}/{chr}.paf.gz"),
		gzi=outpath("minimap2/{chr}/{chr}.paf.gz.gzi")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt *  config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/minimap2:2.31'
	conda:
		'../envs/minimap2.yaml'
	benchmark:
		'benchmarks/{chr}.merge_paf_per_region.benchmark.txt'
	shell:
		'''
		cat {input} > {output.paf}
		bgzip -r {output.paf}
		'''
