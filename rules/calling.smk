from scipy import stats


rule accessible_filter:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_accessible_filter.bed'
    input:
        accessible=rules.accessible_merge.output,
        largeinsert=rules.largeinsert_merge.output,
    params:
        awk=os.path.dirname(workflow.snakefile) + '/tools/accessible_filter.awk'
    shell:
        'bedtools intersect -c -a {input.accessible} -b {input.largeinsert}'
        ' | awk -f {params.awk} > {output}'


rule largeinsert_filter:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_largeinsert_filter.bed'
    input:
        accessible=rules.accessible_filter.output,
        largeinsert=rules.largeinsert_merge.output,
        chrom_size=rules.chrom_sizes.output
    params:
        awk=os.path.dirname(workflow.snakefile) + '/tools/largeinsert_extractor.awk'
    shell:
        'bedtools slop -b 1500 -g {input.chrom_size} -i {input.largeinsert}'
        ' | bedtools intersect -wb -a {input.accessible} -b stdin'
        ' | awk -f {params.awk} > {output}'


def get_depth(wildcards):
    l = float(get_accessible_lambda(wildcards))
    return int(stats.poisson.isf(0.05, l))


rule calling:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_ecDNA.bed'
    input:
        peaks=rules.largeinsert_filter.output,
        bam=rules.merge.output,
        index=rules.index.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/circlehunter.py',
        depth=get_depth
    shell:
        'python {params.script} -d {params.depth} {input.peaks} {input.bam} {output}'