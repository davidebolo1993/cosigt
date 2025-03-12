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

rule add_target_to_queries:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		queries_fasta=config['assemblies'],
		target_fasta=rules.pansnspec_target.output
	output:
		 config['output'] + '/wfmash/queries.fa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/add_target_to_queries.benchmark.txt'
	shell:
		'''
		zcat --force {input.queries_fasta} {input.target_fasta} | \
		awk '/^>/{{f=!d[$1];d[$1]=1}}f' \
		> {output}
		'''

rule samtools_faidx_queries:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.add_target_to_queries.output
	output:
		config['output'] + '/wfmash/queries.fa.fai'
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
		'benchmarks/samtools_faidx_queries.benchmark.txt'
	shell:
		'''
		samtools faidx {input.queries}
		'''                  

rule extract_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		fasta=rules.add_target_to_queries.output,
		fai=rules.samtools_faidx_queries.output
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
		batch='{batch}'
	shell:
		'''
		grep -w {params.batch} {input.fai} | while read f; do
			samtools faidx {input.fasta} $f >> {output}
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
		samtools index {input}
		'''    

rule wfmash_align_batches:
	'''
	https://github.com/waveygang/wfmash
	'''
	input:
		queries_fasta=rules.extract_batches.output
		queries_fai=rules.index_batches.output
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
		'benchmarks/{batch}.merge_batches.benchmark.txt'
	shell:
		'''
		for f in {input}; do
			cat $f >> {output}
		done
		'''