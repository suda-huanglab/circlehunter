include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/accessible.smk'
include: 'rules/largeinsert.smk'
include: 'rules/calling.smk'
include: 'rules/annotate.smk'
include: 'rules/cleanup.smk'
include: 'rules/simulation.smk'


rule all:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/calling/{gsm}_ecDNA_genes.bed'
            for gsm in config['samples']
        ],
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/cleanup.list'
            for gsm in config['samples']
        ]
