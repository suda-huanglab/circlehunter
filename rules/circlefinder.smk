import os


rule calling:
    output:
        bam=config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_microDNA.bam',
        bed=config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_microDNA.bed'
    input:
        bam=rules.merge.output,
        index=rules.index.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/circlehunter.py'
    shell:
        'python {params.script} {input.bam} {output.bed} {output.bam}'
