rule largeinsert_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_tag.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2bed.py',
        mapq=config['params']['mapq'],
        include=config['params']['include'],
        exclude=config['params']['exclude'],
        mismatch=config['params']['mismatch']
    shell:
        'python {params.script} -q {params.mapq} -f {params.include} -F {params.exclude} -r {params.mismatch} -L {input.bed} -i 1500 {input.bam} {output}'


rule largeinsert_pileup:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_pileup.bdg')
    input:
        rules.largeinsert_tag.output
    shell:
        'macs2 pileup --extsize 750 -f BED -i {input} -o {output}'


rule accessible_pileup:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_accessible_pileup.bdg')
    input:
        rules.accessible_tag.output
    shell:
        'macs2 pileup --extsize 750 -f BED -i {input} -o {output}'


rule largeinsert_ratio_value:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_ratio.val')
    input:
        pileup=rules.largeinsert_pileup.output,
        base=rules.accessible_pileup.output
    params:
        awk=os.path.dirname(workflow.snakefile) + '/tools/largeinsert_ratio.awk'
    shell:
        'bedtools intersect -wo -a {input.base} -b {input.pileup} | awk -f {params.awk} > {output}'


def get_largeinsert_ratio(wildcards):
    filename = f'/samples/{wildcards.prefix}/{wildcards.gsm}/largeinsert/{wildcards.gsm}_largeinsert_ratio.val'
    filename = config['workspace'] + filename
    with open(filename) as f:
        return f.read().strip()


rule largeinsert_ratio:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_ratio.bdg')
    input:
        base=rules.accessible_peak.output.lambda_bdg,
        ratio=rules.largeinsert_ratio_value.output
    params:
        ratio=get_largeinsert_ratio
    shell:
        'macs2 bdgopt -i {input.base} -m multiply -p {params.ratio} -o {output}'


def get_accessible_lambda(wildcards):
    filename = f'/samples/{wildcards.prefix}/{wildcards.gsm}/accessible/{wildcards.gsm}_accessible_control_lambda.bdg'
    filename = config['workspace'] + filename
    with open(filename) as f:
        return f.readline().strip().split('\t')[-1]


rule largeinsert_lambda:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_lambda.bdg')
    input:
        largeinsert=rules.largeinsert_ratio.output,
        accessible=rules.accessible_peak.output.lambda_bdg
    params:
        minimum=get_accessible_lambda
    shell:
        'macs2 bdgopt -i {input.largeinsert} -m max -p {params.minimum} -o {output}'


rule largeinsert_pvalue:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_pvalue.bdg')
    input:
        pileup_bdg=rules.largeinsert_pileup.output,
        lambda_bdg=rules.largeinsert_lambda.output
    shell:
        'macs2 bdgcmp -t {input.pileup_bdg} -c {input.lambda_bdg} -m ppois -o {output}'


rule largeinsert_narrowPeak:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_peaks.narrowPeak'
    input:
        rules.largeinsert_pvalue.output
    shell:
        'macs2 bdgpeakcall -i {input} -c 1.301 -l 50 -g 200 -o {output}'


rule largeinsert_merge:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_peaks_merge.bed'
    input:
        peaks=rules.largeinsert_narrowPeak.output,
        chrom_size=rules.chrom_sizes.output
    shell:
        'bedtools sort -g {input.chrom_size} -i {input.peaks}'
        ' | bedtools merge -d 1500 -i stdin -c 4,5 -o first,max > {output}'
 