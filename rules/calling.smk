from scipy import stats


def get_depth(wildcards):
    l = float(get_accessible_lambda(wildcards))
    return int(stats.poisson.isf(0.05, l))


rule calling:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_eccDNA.bed'
    input:
        peaks=rules.largeinsert_accessible.output,
        bam=rules.merge.output,
        index=rules.index.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/circlehunter.py',
        depth=get_depth
    shell:
        'python {params.script} -d {params.depth} {input.peaks} {input.bam} {output}'