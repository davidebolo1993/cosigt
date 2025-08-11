rule ropebwt3_index:
	'''
	https://github.com/lh3/ropebwt3
	- Build index
	'''
	input:
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai
	output:
		fmd=config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fmd',
		ssa=config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fmd.ssa',
		len=config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fmd.len.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['ropebwt3']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['ropebwt3']['time']
	container:
		'docker://davidebolo1993/ropebwt3:3.9'
	conda:
		'../envs/ropebwt3.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.ropebwt3_index.benchmark.txt'
	shell:
		'''
		ropebwt3 build {input.fasta} -t 1 -d -o {output.fmd}
		ropebwt3 ssa {output.fmd} -o {output.ssa} -s 8 
		cut -f 1,2 {input.fai} | gzip > {output.len}
		'''

rule ropebwt3_mem:
	'''
	https://github.com/lh3/ropebwt3
	- Find mem w/ respect to the contigs
	- Convert to .paf
	'''
	input:
		ref_fasta=rules.bedtools_getfasta.output.fasta,
		ref_fai=rules.bedtools_getfasta.output.fai,
		ref_fmd=rules.ropebwt3_index.output.fmd,
		sample_fasta=rules.samtools_fasta_mapped.output
	output:
		config['output'] + '/ropebwt3/{sample}/{chr}/{region}/{region}.realigned.paf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['ropebwt3']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['ropebwt3']['time']
	container:
		'docker://davidebolo1993/ropebwt3:3.9'
	conda:
		'../envs/bwa-mem2.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.ropebwt3_mem.benchmark.txt'
	shell:
		'''
		ropebwt3 mem \
			-l 17 \
			{inpput.fmd} \
			{input.sample_fasta} \
			-p 1 | sh workflow/scripts/convert_mem_to_paf.sh <(cut -f 1,2 {input.fai}) > {output}
		'''