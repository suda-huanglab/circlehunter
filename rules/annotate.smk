import os


rule annotate:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/calling/{gsm}_ecDNA_genes.bed'
    input:
        bed=rules.calling.output,
        db=config['genome']['refgene']
    params:
          script=lambda wildcards: BASE_DIR + '/tools/annotate.py'
    shell:
         'python {params.script} {input.db} {input.bed} {output}'
