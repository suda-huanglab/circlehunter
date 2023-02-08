import os


rule prepare:
    output:
        bam=temp(config['workspace'] + '/samples/{prefix}/{gsm}/circlefinder/{gsm}-hgxx.sorted.bam'),
        index=temp(config['workspace'] + '/samples/{prefix}/{gsm}/circlefinder/{gsm}-hgxx.sorted.bam.bai')
    input:
        bam=rules.merge.output,
        index=rules.index.output
    priority: 20
    shell:
        'ln -s {input.bam} {output.bam} && ln -s {input.index} {output.index}'


rule circlefinder:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/circlefinder/{gsm}-hgxx.microDNA-JT.txt'
    input:
        bam=rules.prepare.output.bam,
        index=rules.prepare.output.index
    priority: 30
    log:
        config['workspace'] + '/samples/{prefix}/{gsm}/log/{gsm}_circlefinder.log'
    benchmark:
        config['workspace'] + '/samples/{prefix}/{gsm}/benchmark/{gsm}_circlefinder.txt'
    params:
        script=lambda wildcards: BASE_DIR + '/tools/circlefinder-pipeline-bwa-mem-samblaster.sh',
        workdir=config['workspace'] + '/samples/{prefix}/{gsm}/circlefinder'
    shell:
        'cd {params.workdir}'
        ' && {params.script} 1 hgxx.fa XX_1.fq.gz XX_2.fq.gz 10 {wildcards.gsm} hgxx > {log} 2>&1'


rule generate_circlefinder_results:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/circlefinder/{gsm}-hgxx.microDNA-JT.txt'
            for gsm in config['samples']
        ]
