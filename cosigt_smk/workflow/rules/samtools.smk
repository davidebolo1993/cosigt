def find_fasta_files(wildcards):
	import glob
	fasta_extensions = ['fa', 'fasta', 'fasta.gz', 'fa.gz']
	base_path = f'resources/alleles/{wildcards.chr}/{wildcards.region}/{wildcards.region}'
	for ext in fasta_extensions:
		pattern = f'{base_path}.{ext}'
		found_files = glob.glob(pattern)
		if found_files:
			return found_files[0]

rule copy_fasta_over:
	'''
	https://github.com/davidebolo1993/cosigt
	- This is necessary to copy the alleles to a standard location
	- Check their extensions
	- If alleles are .gz, assume they are bgzip-compressed - since this is checked by organize.py
	- And potentially compress
	'''
	input:
		find_fasta_files
	output:
		fasta=config['output'] + '/alleles/{chr}/{region}/{region}.fasta.gz',
		fai=config['output'] + '/alleles/{chr}/{region}/{region}.fasta.gz.fai',
		gzi=config['output'] + '/alleles/{chr}/{region}/{region}.fasta.gz.gzi'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['mid']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
	conda:
		'../envs/samtools.yaml'	
	benchmark:
		'benchmarks/{chr}.{region}.copy_fasta_over.benchmark.txt'
	params:
		outdir=config['output'] + '/alleles/{chr}/{region}/'
	shell:
		'''
		if [[ {input} =~ \.gz$ ]]; then
			cp {input} {output.fasta}
		else
			bgzip -c {input} > {output.fasta}
		fi
		samtools faidx {output.fasta}
		'''

rule samtools_fasta_mapped:
	'''
	https://github.com/samtools/samtools
	- Extract reads aligned to the reference
	- Sort them by name (for smart pairing)
	- Convert to fasta
	- Compress
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_alignment_bed.output,
		fasta=config['reference']
	output:
		config['output'] + '/samtools/fasta/{sample}/{chr}/{region}/{region}.mapped.fasta.gz'
	threads:
		config['samtools']['fasta_mapped']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['fasta_mapped']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
	conda:
		'../envs/samtools.yaml'	
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.samtools_fasta_mapped.benchmark.txt'
	params:
		tmpfile=config['output'] + '/samtools/fasta/{sample}/{chr}/{region}/{region}'
	shell:
		'''
		samtools view \
			-T {input.fasta} \
			-@ {threads} \
			-L {input.bed} \
			-M \
			-b \
			{input.sample} | \
			samtools sort \
			-n \
			-@ {threads} \
			-T {params.tmpfile} \
			- | \
			samtools fasta \
			-@ {threads} \
			- | gzip > {output}
		'''

rule samtools_fasta_unmapped:
	'''
	https://github.com/samtools/samtools
	- Extract unmapped reads from each sample alignment
	- Generates an interleaved .fasta
	- Compress
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		fasta=config['reference']
	output:
		config['output'] + '/samtools/fasta/{sample}/unmapped.fasta.gz'
	threads:
		config['samtools']['fasta_unmapped']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['fasta_unmapped']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['fasta_unmapped']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
	conda:
		'../envs/samtools.yaml'	
	benchmark:
		'benchmarks/{sample}.samtools_fasta_unmapped.benchmark.txt'
	params:
		tmpfile=config['output'] + '/samtools/fasta/{sample}/unmapped'
	shell:
		'''
		samtools view \
			-u \
			-f 4 \
			-@ {threads} \
			-T {input.fasta} \
			{input.sample} | \
			samtools sort \
			-n \
			-@ {threads} \
			-T {params.tmpfile} | \
			samtools fasta \
			-0 /dev/null \
			-@ {threads} \
			- | \
			gzip > {output}
		'''
