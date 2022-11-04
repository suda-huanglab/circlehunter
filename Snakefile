include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/accessible.smk'
include: 'rules/largeinsert.smk'
include: 'rules/calling.smk'
include: 'rules/annotate.smk'
include: 'rules/cleanup.smk'
include: 'rules/circlefinder.smk'

if 'simulation' in config:
    include: 'rules/simulation.smk'


import os

BASE_DIR = os.path.dirname(workflow.snakefile)


rule all:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/calling/{gsm}_ecDNA_genes.bed'
            for gsm in config['samples']
        ]
