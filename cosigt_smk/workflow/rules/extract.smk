from glob import glob

rule samtools_view:
	'''
	Samtools view to extract the region
	'''
	input:
		lambda wildcards: glob('resources/cram/{sample}.*am'.format(sample=wildcards.sample))
	output:
		"results/cosigt_results/{sample}/{sample}.region.bam"
	threads:
		5
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		ref=config['reference'],
		region=config['region']
	shell:
		'''
		samtools view -O bam \
		-o {output} \
		-T {params.ref} \
		-@ {threads} \
		{input} \
		{params.region}
		'''

rule samtools_fastq:
	'''
	Samtools fastq to get fastq files from BAM. Ignored in the new version because we are using fasta
	'''
	input:
		rules.samtools_view.output
	output:
		"results/cosigt_results/{sample}/{sample}.region.fastq.gz"
	threads:
		5
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		samtools fastq \
		-@ {threads} \
		{input} | pigz > {output}
		'''

rule samtools_fasta:
        '''
        Samtools fasta to get fasta files from BAM
        '''
        input:
                rules.samtools_view.output
        output:
                "results/cosigt_results/{sample}/{sample}.region.fasta.gz"
        threads:
                5
        container:
                'docker://davidebolo1993/graph_genotyper:latest'
        shell:
                '''
                samtools fasta \
                -@ {threads} \
                {input} | pigz > {output}
                '''
