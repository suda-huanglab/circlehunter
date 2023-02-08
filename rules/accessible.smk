rule accessible_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_tag.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    priority: 20
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_accessible_tag.txt'
    params:
        script=lambda wildcards: BASE_DIR + '/tools/bam2bed.py',
        mapq=config['params']['mapq'],
        include=config['params']['include'],
        exclude=config['params']['exclude'],
        mismatch=config['params']['mismatch']
    shell:
        'python {params.script} -q {params.mapq} -f {params.include} -F {params.exclude} -r {params.mismatch} -L {input.bed} {input.bam} {output}'


rule accessible_peak:
    output:
        peak=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks.narrowPeak',
        treat_bdg=temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_treat_pileup.bdg'),
        lambda_bdg=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_control_lambda.bdg'
    input:
        rules.accessible_tag.output
    priority: 30
    log:
        config['workspace'] + '/samples/{prefix}/{gsm}/log/{gsm}_accessible_peak.log'
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_accessible_peak.txt'
    params:
        genome_size=config['genome']['size'],
        outdir=config['workspace'] + '/samples/{prefix}/{gsm}/accessible',
        name='{gsm}_accessible'
    shell:
        'macs2 callpeak -t {input} -f BEDPE -g {params.genome_size} --keep-dup all --outdir {params.outdir} -n {params.name}'
        ' -B --nomodel -p 0.05 --nolambda --max-gap 200 --min-length 1000 2> {log}'

rule accessible_merge:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks_merge.bed'
    input:
        accessible=rules.accessible_peak.output.peak,
        chrom_size=rules.chrom_sizes.output
    priority: 40
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_accessible_merge.txt'
    shell:
        'bedtools sort -g {input.chrom_size} -i {input.accessible}'
        ' | bedtools merge -d 12500 -i stdin -c 4,5,7 -o first,max,mean > {output}'


rule cutsites:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_cutsites.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=lambda wildcards: BASE_DIR + '/tools/cutsites.py',
        include=lambda wildcards: config['params']['include'],
        exclude=lambda wildcards: config['params']['exclude'],
        mapq=lambda wildcards: config['params']['mapq']
    shell:
        'python {params.script} -f {params.include} -F {params.exclude} -q {params.mapq}'
        ' -L {input.bed} {input.bam} {output}'


rule sort_cutsites:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_cutsites.sorted.bed')
    input:
        rules.cutsites.output
    shell:
        'sort -k1,1V -k2,2n -k3,3n {input} > {output}'


rule bgzip_cutsites:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_cutsites.sorted.bed.gz'
    input:
        rules.sort_cutsites.output
    shell:
        'bgzip -c {input} > {output}'


rule cutsites_index:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_cutsites.sorted.bed.gz.tbi'
    input:
        rules.bgzip_cutsites.output
    shell:
        'tabix {input}'


rule generate_cutsites:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/accessible/{gsm}_cutsites.sorted.bed.gz'
            for gsm in config['samples']
        ] + [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/accessible/{gsm}_cutsites.sorted.bed.gz.tbi'
            for gsm in config['samples']
        ]
