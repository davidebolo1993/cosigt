rule plot_benchmark_time_mem:
    '''
    https://github.com/davidebolo1993/cosigt
    '''
    input:
        expand(config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv', sample=df['sample_id'].tolist(), region=config['region'])
    output:
        time=config['output'] + '/benchmark/resources.time.pdf',
        mem=config['output'] + '/benchmark/resources.mem.pdf'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
        time=lambda wildcards, attempt: attempt * config['default']['time']
    container:
        'docker://davidebolo1993/renv:4.3.3'
    conda:
        '../envs/r.yaml'
    params:
        benchmark_dir='benchmarks'
    shell:
        '''
        Rscript \
        workflow/scripts/plotbenchmark.r \
        {params.benchmark_dir} \
        {output.time} \
        {output.mem}
        '''

rule plot_benchmark_tpr_mask:
    '''
    https://github.com/davidebolo1993/cosigt
    '''
    input:
        expand(config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv', sample=df['sample_id'].tolist(), region=config['region'])
    output:
        config['output'] + '/benchmark/tpr.mask.pdf'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
        time=lambda wildcards, attempt: attempt * config['default']['time']
    container:
        'docker://davidebolo1993/renv:4.3.3'
    conda:
        '../envs/r.yaml'
    params:
        benchmark_dir=config['output'] + '/cosigt',
        cluster_dir=config['output'] + '/cluster',
        dissimilarity_dir=config['output'] + '/odgi/dissimilarity'
    shell:
        '''
        Rscript \
        workflow/scripts/plottpr.mask.r \
        {params.benchmark_dir} \
        {params.cluster_dir} \
        {params.dissimilarity_dir} \
        {output}
        '''

rule plot_benchmark_tpr:
    '''
    https://github.com/davidebolo1993/cosigt
    '''
    input:
        expand(config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv', sample=df['sample_id'].tolist(), region=config['region'])
    output:
        config['output'] + '/benchmark/tpr.pdf'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
        time=lambda wildcards, attempt: attempt * config['default']['time']
    container:
        'docker://davidebolo1993/renv:4.3.3'
    conda:
        '../envs/r.yaml'
    params:
        benchmark_dir=config['output'] + '/cosigt',
        cluster_dir=config['output'] + '/cluster',
        dissimilarity_dir=config['output'] + '/odgi/dissimilarity'
    shell:
        '''
        Rscript \
        workflow/scripts/plottpr.r \
        {params.benchmark_dir} \
        {params.cluster_dir} \
        {params.dissimilarity_dir} \
        {output}
        '''