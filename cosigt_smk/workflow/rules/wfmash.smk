rule pansnspec_target:
	'''
	https://github.com/samtools/samtools
	- Extract the reference chromosome from the reference file
	- Convert chromosome name adding PanSN naming
	- Compress
	- Index
	'''
	input:
		config['reference']
	output:
		fasta=config['output'] + '/wfmash/{chr}/{chr}.fasta.gz',
		fai=config['output'] + '/wfmash/{chr}/{chr}.fasta.gz.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
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

checkpoint generate_batches:
	'''
	https://github.com/davidebolo1993/cosigt
	- Generate batches for parallel alignment
	- With assemblies following PanSN spec, each individual is aligned independently
	'''
	input:
		lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
	output:
		directory(config['output'] + '/wfmash/{chr}/batches/ids')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.generate_batches.benchmark.txt'
	shell:
		'''
		bash workflow/scripts/make_wfmash_batches.sh {input} {output}
		'''

def get_batches(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- Recover names to use as widlcards downstream
	'''
	chr=wildcards.chr
	checkpoint_output = checkpoints.generate_batches.get(chr=chr).output[0]
	batch_files = glob(checkpoint_output + '/*txt')
	return [os.path.basename(f).split('.')[0] for f in batch_files]

rule samtools_faidx_batches:
	'''
	https://github.com/samtools/samtools
	- Extract individual contigs from the original assemblies
	- Compress
	- Index
	'''
	input:
		fai=lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
		ids=config['output'] + '/wfmash/{chr}/batches/ids/{batch}.txt'
	output:
		fasta=temp(config['output'] + '/wfmash/{chr}/batches/fasta/{batch}.fasta.gz'),
		fai=temp(config['output'] + '/wfmash/{chr}/batches/fasta/{batch}.fasta.gz.fai'),
		gzi=temp(config['output'] + '/wfmash/{chr}/batches/fasta/{batch}.fasta.gz.gzi')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * 2 * config['default_high']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{batch}.samtools_faidx_batches.benchmark.txt'
	shell:
		'''
		fai=$(echo {input.fai})
		fasta=$(echo "${{fai%.*}}")
		if [ -f {output.fasta} ]; then
			rm {output.fasta}
		fi
		while read f; do
			samtools faidx $fasta $f | bgzip -c >> {output.fasta}
		done < {input.ids}
		samtools faidx {output.fasta}
		'''

rule wfmash_align_batches:
	'''
	https://github.com/waveygang/wfmash
	- Align individual queries (assemblies) to the target (reference chromosome)
	- Compress
	'''
	input:
		target_fasta=rules.pansnspec_target.output.fasta,
		target_fai=rules.pansnspec_target.output.fai,
		queries_fasta=rules.samtools_faidx_batches.output.fasta,
		queries_fai=rules.samtools_faidx_batches.output.fai,
		queries_gzi=rules.samtools_faidx_batches.output.gzi
	output:
		temp(config['output'] + '/wfmash/{chr}/batches/paf/{batch}.paf.gz')
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
			{params.flags}  | bgzip -c > {output}
		rm -rf {params.tmpdir}
		'''

def get_paf_files(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- Delay resolution
	'''
	batches = get_batches(wildcards)
	return expand(
		config['output'] + '/wfmash/{chr}/batches/paf/{batch}.paf.gz',
		chr=wildcards.chr,
		batch=batches
	)

checkpoint merge_paf_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	- Concatenate the paf files for each chromosome together
	- Ensure temp files are cleaned up after merging
	- Index paf since this is required by impg
	'''
	input:
		get_paf_files
	output:
		paf=config['output'] + '/wfmash/{chr}/{chr}.paf.gz',
		gzi=config['output'] + '/wfmash/{chr}/{chr}.paf.gz.gzi'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/wfmash:0.14.0'
	conda:
		'../envs/wfmash.yaml'
	benchmark:
		'benchmarks/{chr}.merge_paf_per_region.txt'
	params:
		batches_tmp=config['output'] + '/wfmash/{chr}/batches/ids'
	shell:
		'''
		cat {input} > {output.paf}
		bgzip -r {output.paf}
		rm -rf {params.batches_tmp}
		'''

def get_merged_paf(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- For some reason (I don't fully understand, honestly) we need to re-evaluate this
	'''
	checkpoint_output = checkpoints.merge_paf_per_region.get(chr=wildcards.chr).output[0]
	return checkpoint_output