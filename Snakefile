include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/calling.smk'
include: 'rules/annotate.smk'


rule all:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/calling/{gsm}_microDNA.annotated.bed'
            for gsm in config['samples']
        ]
