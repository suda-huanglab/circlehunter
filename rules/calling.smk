from scipy import stats


rule accessible_filter:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_accessible_filter.bed'
    input:
        accessible=rules.accessible_merge.output,
        largeinsert=rules.largeinsert_merge.output
    priority: 80
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_accessible_filter.txt'
    params:
        awk=lambda wildcards: BASE_DIR + '/tools/accessible_filter.awk'
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
    priority: 90
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_largeinsert_filter.txt'
    params:
        awk=lambda wildcards: BASE_DIR + '/tools/largeinsert_extractor.awk'
    shell:
        'bedtools slop -b 750 -g {input.chrom_size} -i {input.largeinsert}'
        ' | bedtools intersect -wb -a {input.accessible} -b stdin'
        ' | awk -f {params.awk} > {output}'


def get_depth(wildcards):
    if 'depth' in config['params']:
        return config['params']['depth']
    l = float(get_accessible_lambda(wildcards))
    return int(stats.poisson.isf(0.05, l))


rule calling:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_ecDNA.bed'
    input:
        peaks=rules.largeinsert_filter.output,
        bam=rules.merge.output,
        index=rules.index.output
    priority: 100
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_calling.txt'
    params:
        script=lambda wildcards: BASE_DIR + '/tools/circlehunter.py',
        mapq=config['params']['mapq'],
        include=config['params']['include'],
        exclude=config['params']['exclude'],
        mismatch=config['params']['mismatch'],
        limit=config['params'].get('limit', 100),
    run:
        depth = get_depth(wildcards)
        shell((
            f'python {params.script} -q {params.mapq} -f {params.include}'
            f' -F {params.exclude} -r {params.mismatch} -d {depth}'
            f' -m {params.limit} {input.peaks} {input.bam} {output}'
        ))
