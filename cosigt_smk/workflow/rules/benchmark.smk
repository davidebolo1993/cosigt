rule plot_benchmark_time_mem:
    input:
        expand(config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv', sample=df['sample_id'].tolist(), region=config['region'])
    output:
        time='benchmark/resources.time.pdf',
        mem='benchmark/resources.mem.pdf'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
        time=lambda wildcards, attempt: attempt * config['default']['time']
    container:
        'docker://davidebolo1993/cosigt_workflow:latest'
    conda:
        '../envs/plot.yaml'
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

rule plot_benchmark_tpr:
    input:
        expand(config['output'] + '/cosigt/{sample}/{region}/cosigt_genotype.tsv', sample=df['sample_id'].tolist(), region=config['region'])
    output:
        'benchmark/tpr.pdf'
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
        time=lambda wildcards, attempt: attempt * config['default']['time']
    container:
        'docker://davidebolo1993/cosigt_workflow:latest'
    conda:
        '../envs/plot.yaml'
    params:
        benchmark_dir=config['output'] + '/cosigt',
        cluster_dir=config['output'] + '/cluster'
    shell:
        '''
        Rscript \
        workflow/scripts/plottpr.r \
        {params.benchmark_dir} \
        {params.cluster_dir} \
        {output}
        '''