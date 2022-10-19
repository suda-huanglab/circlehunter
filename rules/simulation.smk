# ecDNA
SEED = 1665331200
NUM = 5
# reads
READ_LENGTH_RANGES = [35, 50, 75, 100]
READ_DEPTH_RANGES = [10, 20, 30, 40, 50]
FRAGMENT_SIZE_MEAN = 250
FRAGMENT_SIZE_STD = 180


checkpoint mock_ecDNA:
    output:
        bed=config['workspace'] + '/simulation/ecDNA/mock_ecDNA.bed',
        fasta=directory(config['workspace'] + '/simulation/ecDNA/mock_ecDNA.fa')
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
        output=config['workspace'] + '/simulation/ecDNA/mock_ecDNA'
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
        config['workspace'] + '/simulation/ecDNA/mock_ecDNA.fa/ecDNA_{no}.fa'
    input:
        config['workspace'] + '/simulation/ecDNA/mock_ecDNA.fa/ecDNA_{no}.bed'
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
        config['workspace'] + '/simulation/ecDNA/mock_ecDNA_merged.fa'
    input:
        get_all_ecDNA_fasta
    shell:
        'cat {input} > {output}'


rule trim_reads:
    output:
        config['workspace'] + '/simulation/chrom_fastq/{prefix}/{gsm}/L{length}/{srr}_L{length}_{r}.fastq.gz'
    input:
        lambda wildcards: config['samples'][wildcards.gsm]['fastq'][wildcards.srr][f'fq{wildcards.r}']
    shell:
        'zcat {input}'
        ' | awk \'{{if(NR % 4 % 2 == 0){{print(substr($0, 0, {wildcards.length}))}}else{{print($0)}}}}\''
        ' | gzip -c > {output}'


def get_all_samples_fastq(wildcards):
    fastq_files = sum([
        [
            (
                config['workspace'] + f'/simulation/chrom_fastq/{gsm[:6]}/'
                f'{gsm}/L{length}/{srr}_L{length}_{r}.fastq.gz'
            )
            for srr in config['samples'][gsm]['fastq']
            for length in READ_LENGTH_RANGES
            for r in (1, 2)
        ]
        for gsm in config['samples']
    ], [])
    return fastq_files


rule mock_reads:
    output:
        temp(config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_1.fq'),
        temp(config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_2.fq')
    input:
        config['workspace'] + '/simulation/ecDNA/mock_ecDNA.fa/ecDNA_{no}.fa'
    params:
        seed=SEED,
        mean=FRAGMENT_SIZE_MEAN,
        std=FRAGMENT_SIZE_STD,
        output=config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_'
    shell:
        'art_illumina -rs {params.seed} -m {params.mean} -s {params.std} -p -na'
        ' -f {wildcards.depth} -l {wildcards.length} -i {input} -o {params.output}'


rule compress_mock_reads:
    output:
        config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_{r}.fq.gz'
    input:
        config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_{r}.fq'
    shell:
        'gzip -c {input} > {output}'


def get_all_ecDNA_fastq(wildcards):
    out = checkpoints.mock_ecDNA.get(**wildcards).output['fasta']
    ecDNAs = glob_wildcards(f'{out}/ecDNA_{{no}}.bed').no
    fastq_files = [
        (
            config['workspace'] + f'/simulation/ecDNA_fastq/'
            f'D{depth}/L{length}/ecDNA_{no}_{r}.fq.gz'
        )
        for depth in READ_DEPTH_RANGES
        for length in READ_LENGTH_RANGES
        for no in ecDNAs
        for r in (1, 2)
    ]
    return fastq_files


use rule trim as trim_chrom_reads with:
    output:
        fq1=temp(config['workspace'] + '/simulation/chrom_preprocess/{prefix}/{gsm}/L{length}/{srr}_L{length}_r1.trimed.fq.gz'),
        fq2=temp(config['workspace'] + '/simulation/chrom_preprocess/{prefix}/{gsm}/L{length}/{srr}_L{length}_r2.trimed.fq.gz'),
        json=config['workspace'] + '/simulation/chrom_qc/{prefix}/{gsm}/L{length}/{srr}_L{length}_fastp.json',
        html=config['workspace'] + '/simulation/chrom_qc/{prefix}/{gsm}/L{length}/{srr}_L{length}_fastp.html'
    input:
        adapter=config['adapter'],
        fq1=config['workspace'] + '/simulation/chrom_fastq/{prefix}/{gsm}/L{length}/{srr}_L{length}_1.fastq.gz',
        fq2=config['workspace'] + '/simulation/chrom_fastq/{prefix}/{gsm}/L{length}/{srr}_L{length}_2.fastq.gz',
    log:
        config['workspace'] + '/simulation/chrom_log/{prefix}/{gsm}/L{length}/{srr}_L{length}_fastp.log'


rule generate_fastq:
    input:
        get_all_samples_fastq,
        get_all_ecDNA_fastq
