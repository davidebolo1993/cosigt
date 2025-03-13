rule pansnspec_target:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		ref=config['reference']
	output:
		config['output'] + '/wfmash/target.fa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	benchmark:
		'benchmarks/pansnspec_toref.benchmark.txt'
	conda:
		'../envs/samtools.yaml'
	params:
		path=config['path']
	shell:
		'''
		samtools faidx \
		{input.ref} \
		$(echo {params.path} | cut -d "#" -f 2) | \
		sed "1 s/^.*$/>{params.path}/" \
		> {output}
		'''

rule samtools_faidx_target:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.pansnspec_target.output
	output:
		config['output'] + '/wfmash/target.fa.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/samtools_faidx_target.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''
   
rule extract_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		config['assemblies']
	output:
		config['output'] + '/wfmash/batches/{batch}.fa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{batch}.extract_batch.benchmark.txt'
	params:
		batch='{batch}',
		fai=config['assemblies'] + '.fai'
	shell:
		'''
		grep -w {params.batch} {params.fai} | cut -f 1 | while read f; do
			samtools faidx {input} $f >> {output}
		done
		'''

rule samtools_faidx_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.extract_batches.output
	output:
		config['output'] + '/wfmash/batches/{batch}.fa.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{batch}.index_batches.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''    

rule wfmash_align_batches:
	'''
	https://github.com/waveygang/wfmash
	'''
	input:
		queries_fasta=rules.extract_batches.output,
		queries_fai=rules.samtools_faidx_batches.output,
		target_fasta=rules.pansnspec_target.output,
		target_fai=rules.samtools_faidx_batches.output
	output:
		config['output'] + '/wfmash/batches/{batch}.paf'
	threads:
		config['wfmash']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['wfmash']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['wfmash']['time']
	container:
		'docker://davidebolo1993/wfmash:0.14.0'
	conda:
		'../envs/wfmash.yaml'
	benchmark:
		'benchmarks/{batch}.wfmash_align_batches.benchmark.txt'
	params:
		flags=config['wfmash']['params'],
		tmpdir=config['wfmash']['tmpdir'] 
	shell:
		'''
		wfmash \
			{input.target_fasta} \
			{input.queries_fasta} \
			-X \
			-t {threads} \
			-B {params.tmpdir} \
			{params.flags} > {output}
		'''

rule merge_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		expand(config['output'] + '/wfmash/batches/{batch}.paf', batch=sorted(batch_set))
	output:
		config['output'] + '/wfmash/queries_to_target.paf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/merge_batches.benchmark.txt'
	shell:
		'''
		for f in {input}; do
			cat $f >> {output}
		done
		'''