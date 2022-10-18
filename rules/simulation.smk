SEED = 1665331200
NUM = 5


checkpoint mock_ecDNA:
    output:
        bed=config['workspace'] + '/simulation/ecDNA/mock-ecDNA.bed',
        fasta=directory(config['workspace'] + '/simulation/ecDNA/mock-ecDNA.fa')
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/mock.py',
        num=NUM,
        multiply=6,
        seed=SEED,
        length_loc=12,
        length_scale=3.5,
        length_minimum=5000,
        segments_loc=1,
        segments_scale=1,
        max_consecutive_n=10000,
        chrom_sizes=config['genome']['chrom_sizes'],
        blacklist=config['genome']['blacklist'],
        fasta=config['genome']['fasta'],
        output=config['workspace'] + '/simulation/ecDNA/mock-ecDNA'
    threads:
        workflow.cores
    benchmark:
        config['workspace'] + '/simulation/benchmarks/mock/mock.tsv'
    shell:
        'python {params.script} -n {params.num} -m {params.multiply}'
        ' -s {params.seed} -p {threads}'
        ' --length-loc {params.length_loc} --length-scale {params.length_scale}'
        ' --length-minimum {params.length_minimum}'
        ' --segments-loc {params.segments_loc} --segments-scale {params.segments_scale}'
        ' --max-consecutive-n {params.max_consecutive_n}'
        ' {params.chrom_sizes} {params.blacklist} {params.fasta} {params.output}'


rule get_fasta:
    output:
        config['workspace'] + '/simulation/ecDNA/mock-ecDNA.fa/ecDNA_{no}.fa'
    input:
        config['workspace'] + '/simulation/ecDNA/mock-ecDNA.fa/ecDNA_{no}.bed'
    params:
        fasta=config['genome']['fasta']
    shell:
        'echo ">ecDNA_{wildcards.no}" > {output}'
        ' && (bedtools getfasta -nameOnly -s -fi {params.fasta} -bed {input}'
        ' | grep -v "^>" | sed ":a; N; s/\\n//; ta" | fold) >> {output}'


def get_all_ecDNA_fasta(wildcards):
    out = checkpoints.mock_ecDNA.get(**wildcards).output['fasta']
    bed_files = expand(
        f'{out}/ecDNA_{{no}}.fa',
        no=glob_wildcards(f'{out}/ecDNA_{{no}}.bed').no
    )
    return bed_files


rule merge_fasta:
    output:
        config['workspace'] + '/simulation/ecDNA/mock-ecDNA-merged.fa'
    input:
        get_all_ecDNA_fasta
    shell:
        'cat {input} > {output}'
