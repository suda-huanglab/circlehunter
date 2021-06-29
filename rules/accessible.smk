rule accessible_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_tag.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2bed.py',
        mapq=config['params']['mapq']
    shell:
        'python {params.script} -q {params.mapq} -L {input.bed} {input.bam} {output}'


rule accessible_peak:
    output:
        peak=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks.narrowPeak',
        treat_bdg=temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_treat_pileup.bdg'),
        lambda_bdg=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_control_lambda.bdg'
    input:
        rules.accessible_tag.output
    params:
        outdir=config['workspace'] + '/samples/{prefix}/{gsm}/accessible',
        name='{gsm}_accessible'
    shell:
        'macs2 callpeak -t {input} -f BEDPE -g hs --keep-dup all --outdir {params.outdir} -n {params.name}'
        ' -B --nomodel -p 0.05 --nolambda --max-gap 200 --min-length 1000'
