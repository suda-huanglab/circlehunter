import os


rule annotate:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_eccDNA_genes.bed'
    input:
        bed=rules.calling.output,
        db=config['genomes'][config['assembly']]['refgene']
    params:
          script=os.path.dirname(workflow.snakefile) + '/tools/annotate.py'
    shell:
         'python {params.script} {input.db} {input.bed} {output}'
