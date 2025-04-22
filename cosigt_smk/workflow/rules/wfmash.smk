rule pansnspec_target:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		config['reference']
	output:
		config['output'] + '/wfmash/{chr}/target.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
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
		sed "1 s/^.*$/>{params.pansn}/" \
		> {output}
		'''

rule samtools_faidx_target:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.pansnspec_target.output
	output:
		config['output'] + '/wfmash/{chr}/target.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.samtools_faidx_target.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''

checkpoint generate_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
	output:
		directory(config['output'] + '/wfmash/{chr}/batches/ids')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	benchmark:
		'benchmarks/{chr}.generate_batches.benchmark.txt'
	shell:
		'''
		bash workflow/scripts/make_wfmash_batches.sh {input} {output}
		'''

def get_batches(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	chr=wildcards.chr
	checkpoint_output = checkpoints.generate_batches.get(chr=chr).output[0]
	batch_files = glob(checkpoint_output + '/*ids')
	return [os.path.basename(f).split('.')[0] for f in batch_files]

rule samtools_faidx_batches:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		fai=lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
		ids=config['output'] + '/wfmash/{chr}/batches/ids/{batch}.ids'
	output:
		config['output'] + '/wfmash/{chr}/batches/fasta/{batch}.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{batch}.samtools_faidx_batches.benchmark.txt'
	shell:
		'''
		fai=$(echo {input.fai})
		fasta=$(echo "${{fai%.*}}")
		touch {output}
		while read f; do
			samtools faidx $fasta $f >> {output}
		done < {input.ids}
		'''

rule samtools_faidx_index_batches:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.samtools_faidx_batches.output
	output:
		config['output'] + '/wfmash/{chr}/batches/fasta/{batch}.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{batch}.samtools_faidx_index_batches.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''

rule wfmash_align_batches:
	'''
	https://github.com/waveygang/wfmash
	'''
	input:
		target_fasta=rules.pansnspec_target.output,
		target_fai=rules.samtools_faidx_target.output,
		queries_fasta=rules.samtools_faidx_batches.output,
		queries_fai=rules.samtools_faidx_index_batches.output
	output:
		config['output'] + '/wfmash/{chr}/batches/paf/{batch}.paf'
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
		'benchmarks/{chr}.{batch}.wfmash_align_batches.benchmark.txt'
	params:
		flags=config['wfmash']['params'],
		tmpdir=config['wfmash']['tmpdir'] + '/{chr}/{batch}'
	shell:
		'''
		mkdir -p {params.tmpdir}
		wfmash \
			{input.target_fasta} \
			{input.queries_fasta} \
			-X \
			-t {threads} \
			-B {params.tmpdir} \
			{params.flags} > {output}
		rm -rf {params.tmpdir}
		'''
