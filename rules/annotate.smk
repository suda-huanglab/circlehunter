import os


rule annotate:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_microDNA.annotated.bed'
    input:
        bed=rules.calling.output.bed,
        db=config['genomes']['hg38']['refgene']
    params:
          script=os.path.dirname(workflow.snakefile) + '/tools/annotate.py'
    shell:
         'python {params.script} {input.db} {input.bed} {output}'
