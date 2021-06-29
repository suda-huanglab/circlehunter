rule largeinsert_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_tag.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2bed.py',
        mapq=config['params']['mapq']
    shell:
        'python {params.script} -q {params.mapq} -L {input.bed} -i 1500 {input.bam} {output}'


rule largeinsert_pileup:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_pileup.bdg')
    input:
        rules.largeinsert_tag.output
    shell:
        'macs2 pileup -i {input} -f BEDPE -o {output}'


rule largeinsert_ratio_value:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_ratio.val')
    input:
        pileup=rules.largeinsert_pileup.output,
        base=rules.accessible_peak.output.lambda_bdg
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
        'macs2 bdgpeakcall -i {input} -c 1.301 -l 50 -g 1500 -o {output}'


rule largeinsert_slop:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_peaks_slop.bed')
    input:
        peaks=rules.largeinsert_narrowPeak.output,
        chrom_size=rules.chrom_sizes.output
    shell:
        'bedtools slop -b 1500 -g {input.chrom_size} -i {input.peaks}'
        ' | bedtools sort -g {input.chrom_size} -i  stdin'
        ' | bedtools merge -i stdin > {output}'


rule largeinsert_accessible:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_accessible.bed'
    input:
        largeinsert=rules.largeinsert_slop.output,
        accessible=rules.accessible_peak.output.peak,
        chrom_size=rules.chrom_sizes.output
    params:
        accessible_filter=os.path.dirname(workflow.snakefile) + '/tools/accessible_filter.awk',
        largeinsert_extractor=os.path.dirname(workflow.snakefile) + '/tools/largeinsert_extractor.awk'
    shell:
        'bedtools sort -g {input.chrom_size} -i {input.accessible}'
        ' | bedtools merge -d 12500 -i stdin -c 4,5 -o first,first '
        ' | bedtools intersect -c -a stdin -b {input.largeinsert}'
        ' | awk -f {params.accessible_filter}'
        ' | bedtools intersect -wb -a stdin -b {input.largeinsert}'
        ' | awk -f {params.largeinsert_extractor} > {output}'
 