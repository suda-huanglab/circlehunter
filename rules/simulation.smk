## get params from config file
# ecDNA
SEED = config['simulation']['seed']
NUM = config['simulation']['num']
REPEAT = config['simulation']['repeat']

# reads
FRAGMENT_SIZE_MEAN = config['simulation']['fragment_size_mean']
FRAGMENT_SIZE_STD = config['simulation']['fragment_size_std']

# depth test
READ_DEPTH_RANGES = config['simulation']['read_depth_ranges']
READ_DEPTH_TEST_LENGTH = config['simulation']['read_depth_test_length']

# length test
READ_LENGTH_RANGES = config['simulation']['read_length_ranges']
READ_LENGTH_TEST_DEPTH = config['simulation']['read_length_test_depth']


# segments
SEGMENTS_LOC = config['simulation'].get('segments_loc', 1)
SEGMENTS_SCALE = config['simulation'].get('segments_scale', 1)


checkpoint mock_ecDNA_regions:
    output:
        bed=config['workspace'] + '/simulation/ecDNA/mock_ecDNA.bed',
        fasta=directory(config['workspace'] + '/simulation/ecDNA/mock_ecDNA.fa')
    params:
        script=lambda wildcards: BASE_DIR + '/tools/mock.py',
        num=NUM * REPEAT,
        multiply=10,
        seed=SEED,
        length_loc=12,
        length_scale=3.5,
        length_minimum=5000,
        length_maximum=500_000,
        segments_loc=SEGMENTS_LOC,
        segments_scale=SEGMENTS_SCALE,
        max_consecutive_n=10000,
        mock_regions=config['genome']['mock_regions'],
        output=config['workspace'] + '/simulation/ecDNA/mock_ecDNA'
    threads:
        workflow.cores
    benchmark:
        config['workspace'] + '/simulation/benchmarks/mock/mock.tsv'
    shell:
        'python {params.script} -n {params.num} -m {params.multiply} -s {params.seed}'
        ' --length-loc {params.length_loc} --length-scale {params.length_scale}'
        ' --length-minimum {params.length_minimum} --length-maximum {params.length_maximum}'
        ' --segments-loc {params.segments_loc} --segments-scale {params.segments_scale}'
        ' {params.mock_regions} {params.output}'


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
    out = checkpoints.mock_ecDNA_regions.get(**wildcards).output['fasta']
    bed_files = expand(
        f'{out}/ecDNA_{{no}}.fa',
        no=range(1, NUM * REPEAT + 1)
    )
    return bed_files


rule mock_ecDNA:
    input:
        get_all_ecDNA_fasta


rule trim_reads:
    output:
        temp(config['workspace'] + '/simulation/chrom_fastq/{prefix}/{gsm}/L{length}/{srr}_L{length}_{r}.fastq.gz')
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
    log:
        config['workspace'] + '/simulation/ecDNA_log/D{depth}/L{length}/ecDNA_{no}_art.log'
    params:
        seed=SEED,
        mean=FRAGMENT_SIZE_MEAN,
        std=FRAGMENT_SIZE_STD,
        output=config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_'
    shell:
        'art_illumina -rs {params.seed} -m {params.mean} -s {params.std} -p -na'
        ' -f {wildcards.depth} -l {wildcards.length} -i {input} -o {params.output} > {log} 2>&1'


rule compress_mock_reads:
    output:
        temp(config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_{r}.fq.gz')
    input:
        config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_{r}.fq'
    shell:
        'gzip -c {input} > {output}'


def get_all_ecDNA_fastq(wildcards):
    out = checkpoints.mock_ecDNA_regions.get(**wildcards).output['fasta']
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


use rule mapping as mapping_chrom_reads with:
    output:
        config['workspace'] + '/simulation/chrom_mapping/{prefix}/{gsm}/L{length}/{srr}_L{length}.sorted.bam'
    input:
        fq1=rules.trim_chrom_reads.output.fq1,
        fq2=rules.trim_chrom_reads.output.fq2
    priority: 100
    log:
       bwa=config['workspace'] + '/simulation/chrom_log/{prefix}/{gsm}/L{length}/{srr}_L{length}_bwa.log',
       samblaster=config['workspace'] + '/simulation/chrom_log/{prefix}/{gsm}/L{length}/{srr}_L{length}_samblaster.log'
    params:
        rg='\'@RG\\tID:{srr}_L{length}\\tSM:{gsm}\\tLB:{srr}_L{length}\\tPL:ILLUMINA\'',
        index=config['genome']['bwa_index'],
        tmp=config['workspace'] + '/simulation/chrom_mapping/{prefix}/{gsm}/L{length}/{srr}_L{length}.tmp'


def get_all_samples_bam(wildcards):
    bam_files = [
        (
            config['workspace'] + '/simulation/chrom_mapping/'
            f'{gsm[:6]}/{gsm}/L{length}/{srr}_L{length}.sorted.bam'
        )
        for gsm in config['samples']
        for srr in config['samples'][gsm]['fastq']
        for length in READ_LENGTH_RANGES
    ]
    return bam_files


use rule trim as trim_ecNDA_reads with:
    output:
        fq1=temp(config['workspace'] + '/simulation/ecDNA_preprocess/D{depth}/L{length}/ecDNA_{no}_1.fq.gz'),
        fq2=temp(config['workspace'] + '/simulation/ecDNA_preprocess/D{depth}/L{length}/ecDNA_{no}_2.fq.gz'),
        json=config['workspace'] + '/simulation/ecDNA_qc/D{depth}/L{length}/ecDNA_{no}_fastp.json',
        html=config['workspace'] + '/simulation/ecDNA_qc/D{depth}/L{length}/ecDNA_{no}_fastp.html'
    input:
        adapter=config['adapter'],
        fq1=config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_1.fq.gz',
        fq2=config['workspace'] + '/simulation/ecDNA_fastq/D{depth}/L{length}/ecDNA_{no}_2.fq.gz',
    log:
        config['workspace'] + '/simulation/ecDNA_log/D{depth}/L{length}/ecDNA_{no}_fastp.log'


use rule mapping as mapping_ecDNA_reads with:
    output:
        config['workspace'] + '/simulation/ecDNA_mapping/D{depth}/L{length}/ecDNA_{no}.sorted.bam'
    input:
        fq1=rules.trim_ecNDA_reads.output.fq1,
        fq2=rules.trim_ecNDA_reads.output.fq2
    priority: 50
    log:
       bwa=config['workspace'] + '/simulation/ecDNA_log/D{depth}/L{length}/ecDNA_{no}_bwa.log',
       samblaster=config['workspace'] + '/simulation/ecDNA_log/D{depth}/L{length}/ecDNA_{no}_samblaster.log'
    params:
        rg='\'@RG\\tID:ecDNA_{no}_D{depth}_L{length}\\tSM:ecDNA_{no}\\tLB:ecDNA_{no}_D{depth}_L{length}\\tPL:ILLUMINA\'',
        index=config['genome']['bwa_index'],
        tmp=config['workspace'] + '/simulation/ecDNA_mapping/D{depth}/L{length}/ecDNA_{no}.tmp'


def get_all_test_samples(wildcards):
    from itertools import groupby, product
    out = checkpoints.mock_ecDNA_regions.get(**wildcards).output['fasta']
    ecDNAs = list(sorted(map(int, glob_wildcards(f'{out}/ecDNA_{{no}}.bed').no)))[:NUM * REPEAT]
    ecDNAs = {group: list(batch) for group, batch in groupby(ecDNAs, lambda x: (x - 1) // NUM)}
    depth_test_samples = {
        f'{gsm}_D{depth}_L{READ_DEPTH_TEST_LENGTH}_ecDNA_G{group}': {
            'bam': {
                **{
                    # chrom bam files
                    f'{srr}_L{READ_DEPTH_TEST_LENGTH}': (
                        config['workspace'] + f'/simulation/chrom_mapping/{gsm[:6]}/{gsm}/'
                        f'L{READ_DEPTH_TEST_LENGTH}/{srr}_L{READ_DEPTH_TEST_LENGTH}.sorted.bam'
                    )
                    for srr in config['samples'][gsm]['fastq']
                }, **{
                    f'ecDNA_{no}': (
                        config['workspace'] + f'/simulation/ecDNA_mapping/D{depth}/'
                        f'L{READ_DEPTH_TEST_LENGTH}/ecDNA_{no}.sorted.bam'
                    )
                    for no in batch
                }
            }
        }
        for gsm, depth, (group, batch) in product(
            config['samples'], READ_DEPTH_RANGES, ecDNAs.items()
        )
    }
    length_test_samples = {
        f'{gsm}_D{READ_LENGTH_TEST_DEPTH}_L{length}_ecDNA_G{group}': {
            'bam': {
                **{
                    # chrom bam files
                    f'{srr}_L{length}': (
                        config['workspace'] + f'/simulation/chrom_mapping/{gsm[:6]}/{gsm}/'
                        f'L{length}/{srr}_L{length}.sorted.bam'
                    )
                    for srr in config['samples'][gsm]['fastq']
                }, **{
                    f'ecDNA_{no}': (
                        config['workspace'] + f'/simulation/ecDNA_mapping/D{READ_LENGTH_TEST_DEPTH}/'
                        f'L{length}/ecDNA_{no}.sorted.bam'
                    )
                    for no in batch
                }
            }
        }
        for gsm, length, (group, batch) in product(
            config['samples'], READ_LENGTH_RANGES, ecDNAs.items()
        )
    }
    samples = {**depth_test_samples, **length_test_samples}
    return samples


def get_all_bam(wildcards):
    samples = get_all_test_samples(wildcards)
    bam_files = sum([
        list(sample['bam'].values()) for sample in samples.values()
    ], [])
    return bam_files


rule simulation_config:
    output:
        config['workspace'] + '/simulation/simulation.yaml'
    input:
        get_all_bam
    run:
        import json, yaml
        simulation = json.loads(json.dumps(config))
        simulation['workspace'] = config['workspace'] + '/simulation'
        simulation['samples'] = get_all_test_samples(wildcards)
        with open(output[0], 'w') as f:
            yaml.dump(simulation, f)
