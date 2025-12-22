rule meryl_build_reference_db:
	'''
	https://github.com/marbl/meryl
	- Count all the 31-mers in the reference
	- Output to a meryl db
	'''
	input:
		config['reference']
	output:
		config['output'] + '/meryl/reference.done'
	threads:
		20
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['meryl']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['meryl']['time']
	container:
		'docker://davidebolo1993/kfilt:0.1.1'
	benchmark:
		'benchmarks/meryl_buil_reference_db.benchmark.txt'
	conda:
		'../envs/kfilt.yaml'
	params:
		outdir=config['output'] + '/meryl/reference'
	shell:
		'''
		if [ -d {params.outdir} ]; then
			rm -rf {params.outdir}
		fi
		memgb=$(echo {resources.mem_mb} | awk '{{print int(1 + $1/1024 - 0.000001)}}')
		meryl count \
			k=31 \
			threads={threads} \
			memory=$memgb \
			{input} \
			output \
			{params.outdir} \
		&& touch {output}
		'''

rule meryl_build_alleles_db_meryl_difference_kfilt_index:
	'''
	https://github.com/marbl/meryl
	- Count all the 31-mers in the alleles
	- Only retain 31-mers in the alleles that are not in the reference
	- Build the index of the 31-mers with kfilt
	'''
	input:
		db=rules.meryl_build_reference_db.output,
		fasta=rules.bedtools_getfasta.output.fasta
	output:
		config['output'] + '/kfilt/index/{chr}/{region}/{region}.kfilt.idx'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['high']['time']
	container:
		'docker://davidebolo1993/kfilt:0.1.1'
	benchmark:
		'benchmarks/{chr}.{region}.meryl_build_alleles_db_meryl_difference_kfilt_index.benchmark.txt'
	conda:
		'../envs/kfilt.yaml'
	params:
		outdir_alleles=config['output'] + '/meryl/{chr}/{region}/{region}',
		outdir_diff=config['output'] + '/meryl/{chr}/{region}/{region}_unique',
		out_diff=config['output'] + '/meryl/{chr}/{region}/{region}.unique_kmers.txt',
		indir=config['output'] + '/meryl/reference'
	shell:
		'''
		if [ -d {params.outdir_alleles} ]; then
			rm -rf {params.outdir_alleles}
		fi
		if [ -d {params.outdir_diff} ]; then
			rm -rf {params.outdir_diff}
		fi
		memgb=$(echo {resources.mem_mb} | awk '{{print int(1 + $1/1024 - 0.000001)}}')
		mkdir -p {params.outdir_alleles}
		meryl count \
			k=31 \
			threads={threads} \
			memory=$memgb \
			{input.fasta} \
			output \
			{params.outdir_alleles}
		mkdir -p {params.outdir_diff}
		meryl difference \
			{params.outdir_alleles} \
			{params.indir} \
			output \
			{params.outdir_diff}
		rm -rf {params.outdir_alleles}
		meryl print {params.outdir_diff} > {params.out_diff}
		rm -rf {params.outdir_diff}
		kfilt \
			build \
			-k {params.out_diff} \
			-K 31 \
			-o {output}
		'''

rule kfilt_filter_unmapped:
	'''
	https://github.com/davidebolo1993/kfilt
	- Filter unmapped reads per sample
	- Keep only those having at least one 31-mer matching the index
	'''
	input:
		idx=rules.meryl_build_alleles_db_meryl_difference_kfilt_index.output,
		sample=rules.samtools_fasta_unmapped.output
	output:
		config['output'] + '/kfilt/{sample}/{chr}/{region}/{region}.unmapped.fasta.gz'
	threads:
		config['kfilt']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['kfilt']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['kfilt']['time']
	container:
		'docker://davidebolo1993/kfilt:0.1.1'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.kfilt_filter_unmapped.benchmark.txt'
	conda:
		'../envs/kfilt.yaml'
	shell:
		'''
		kfilt \
			filter \
			-I {input.sample} \
			-o {output} \
			-f fasta \
			-z \
			-i {input.idx} \
			-n 1 \
			-m 0 \
			-t {threads}
		'''

rule combine_mapped_unmapped:
	'''
	https://github.com/davidebolo1993/cosigt
	- Concatenate unmapped and mapped reads
	'''
	input:
		fasta_mapped=rules.samtools_fasta_mapped.output,
		fasta_unmapped=rules.kfilt_filter_unmapped.output
	output:
		config['output'] + '/combine/{sample}/{chr}/{region}/{region}.fasta.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['mid']['time']
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.combine_mapped_unmapped.benchmark.txt'
	shell:
		'''
		cat \
			{input.fasta_mapped} \
			{input.fasta_unmapped} > {output}
		'''
	
