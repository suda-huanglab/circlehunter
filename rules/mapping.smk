import os


rule mapping:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{srr}.sorted.bam')
    input:
        fq1=rules.trim.output.fq1,
        fq2=rules.trim.output.fq2
    priority: 10
    log:
       bwa=config['workspace'] + '/samples/{prefix}/{gsm}/log/{srr}_bwa.log',
       samblaster=config['workspace'] + '/samples/{prefix}/{gsm}/log/{srr}_samblaster.log'
    params:
        rg='\'@RG\\tID:{srr}\\tSM:{gsm}\\tLB:{srr}\\tPL:ILLUMINA\'',
        index=config['genome']['bwa_index'],
        tmp=config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{srr}.tmp',
    threads: 8 if workflow.cores > 8 else workflow.cores
    shell:
        'bwa mem -t {threads} -R {params.rg} {params.index} {input.fq1} {input.fq2} 2>{log.bwa}'
        ' | samblaster 2>{log.samblaster}'
        ' | samtools sort -T {params.tmp} -o {output} -'


def get_all_bam(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.gsm[:6]}/{wildcards.gsm}/' \
                              f'mapping/{srr}.sorted.bam'
        for srr in config['samples'][wildcards.gsm]
    ]


rule merge:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{gsm}.sorted.bam'
    input:
        get_all_bam
    shell:
        'samtools merge {output} {input}'


rule index:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{gsm}.sorted.bam.bai'
    input:
        rules.merge.output
    priority: 5
    shell:
        'samtools index {input}'


rule chrom_sizes:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{gsm}.chrom_sizes'
    input:
        bam=rules.merge.output,
        index=rules.index.output
    params:
        regex=os.path.dirname(workflow.snakefile) + '/tools/chrom_size.regex',
        awk=os.path.dirname(workflow.snakefile) + '/tools/chrom_size.awk'
    shell:
        'samtools idxstats {input.bam} | grep -f {params.regex} | awk -f {params.awk} > {output}'


rule clean_bed:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/mapping/{gsm}_clean.bed'
    input:
        rules.chrom_sizes.output
    params:
        blacklist=config['genome']['blacklist']
    shell:
        'bedtools complement -g {input} -i {params.blacklist} > {output}'
