rule accessible_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_tag.bed')
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
        'python {params.script} -q {params.mapq} -f {params.include} -F {params.exclude} -r {params.mismatch} -L {input.bed} {input.bam} {output}'


rule accessible_peak:
    output:
        peak=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks.narrowPeak',
        treat_bdg=temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_treat_pileup.bdg'),
        lambda_bdg=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_control_lambda.bdg'
    input:
        rules.accessible_tag.output
    params:
        genome_size=config['genome']['size'],
        outdir=config['workspace'] + '/samples/{prefix}/{gsm}/accessible',
        name='{gsm}_accessible'
    shell:
        'macs2 callpeak -t {input} -f BEDPE -g {params.genome_size} --keep-dup all --outdir {params.outdir} -n {params.name}'
        ' -B --nomodel -p 0.05 --nolambda --max-gap 200 --min-length 1000'

rule accessible_merge:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks_merge.bed'
    input:
        accessible=rules.accessible_peak.output.peak,
        chrom_size=rules.chrom_sizes.output
    shell:
        'bedtools sort -g {input.chrom_size} -i {input.accessible}'
        ' | bedtools merge -d 12500 -i stdin -c 4,5 -o first,max > {output}'
