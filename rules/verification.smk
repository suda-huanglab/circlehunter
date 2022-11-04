from functools import lru_cache
import os


FOLDCHANGE_RANGES = 5, 51, 5
READ_LENGTH_FOLDCHANGES = [30, 50]
READ_LENGTH_RANGES = 35, 101, 5

wildcard_constraints:
    foldchange="\d+"


rule verification_config:
    output:
        config['workspace'] + '/verification/verify.yaml'
    input:
        fasta=config['workspace'] + '/verification/genome/mock.fa',
        blacklist=config['workspace'] + '/verification/genome/mock.blacklist.bed',
        mock=config['workspace'] + '/verification/genome/mock.ecDNA.bed',
        fai=config['workspace'] + '/verification/genome/mock.fa.fai',
        bwt=config['workspace'] + '/verification/genome/mock.fa.bwt',
        reads=[
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{gsm}_chr22.list'
            for gsm in config['samples']
        ],
        fastq=[
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_{chrom}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for chrom in ('chr22', 'chrM')
            for r in (1, 2)
        ],
        depth=[
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{gsm}_chrM_depth.value'
            for gsm in config['samples']
        ],
        all_fastq=[
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_chrM_{foldchange}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for foldchange in range(*FOLDCHANGE_RANGES)
            for r in (1, 2)
        ] + [
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_chr22_L{length}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for foldchange in READ_LENGTH_FOLDCHANGES
            for length in range(*READ_LENGTH_RANGES)
            for r in (1, 2)
        ] + [
            config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_chrM_{foldchange}_L{length}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for foldchange in READ_LENGTH_FOLDCHANGES
            for length in range(*READ_LENGTH_RANGES)
            for r in (1, 2)
        ]
    run:
        from pyfaidx import Fasta
        import pandas as pd
        import json
        import yaml
        genome = Fasta(input.fasta)
        genome = {
            'size': len(genome['chr23']),
            'fasta': input.fasta,
            'bwa_index': input.fasta,
            'blacklist': input.blacklist,
            'refgene': config['genome']['refgene']
        }
        chroms = [f'chrM_{foldchange}' for foldchange in range(*FOLDCHANGE_RANGES)] + ['chrM']
        samples1 = {
            f'{gsm}_{foldchange}': {
                f'{srr}_{chrom}': {
                    f'fq{r}': config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_{chrom}_{r}.fastq.gz'
                    for r in (1, 2)
                }
                for srr in config['samples'][gsm]
                for chrom in ('chr22', f'chrM_{foldchange}')
            }
            for gsm in config['samples']
            for foldchange in range(*FOLDCHANGE_RANGES)
        }
        samples2 = {
            f'{gsm}_{foldchange}_L{length}': {
                f'{srr}_{chrom}_L{length}': {
                    f'fq{r}': config['workspace'] + f'/verification/fastq/{gsm[:6]}/{gsm}/{srr}_{chrom}_L{length}_{r}.fastq.gz'
                    for r in (1, 2)
                }
                for srr in config['samples'][gsm]
                for chrom in ('chr22', f'chrM_{foldchange}')
            }
            for gsm in config['samples']
            for foldchange in READ_LENGTH_FOLDCHANGES
            for length in range(*READ_LENGTH_RANGES)
        }
        samples = {**samples1, **samples2}
        verification = {
            'genome': genome,
            'adapter': json.loads(json.dumps(config['adapter'])),
            'params': json.loads(json.dumps(config['params'])),
            'workspace': config['workspace'] + '/verification',
            'samples': samples
        }
        with open(output[0], 'w') as f:
            yaml.dump(verification, f)


rule mock_genome:
    output:
        fasta=config['workspace'] + '/verification/genome/mock.fa',
        blacklist=config['workspace'] + '/verification/genome/mock.blacklist.bed',
        mock=config['workspace'] + '/verification/genome/mock.ecDNA.bed'
    input:
        fasta=config['genome']['fasta'],
        blacklist=config['genome']['blacklist']
    run:
        from pyfaidx import Fasta
        import pandas as pd
        genome = Fasta(input.fasta)
        l = len(genome['chr22']) // 3
        s = len(genome['chrM']) // 2
        f = open(f'{output.fasta}.tmp', 'w')
        f.write('>chr23\n')
        f.write(genome['chr22'][:l].seq)
        f.write(genome['chrM'][:s].seq)
        f.write(genome['chr22'][l:2*l].seq)
        f.write((-genome['chrM'][s:]).seq)
        f.write(genome['chr22'][2*l:].seq)
        f.close()
        shell(f'fold {output.fasta}.tmp > {output.fasta} && rm {output.fasta}.tmp')
        blacklist = pd.read_table(input.blacklist, names=['chrom', 'start', 'end'])
        blacklist = blacklist[blacklist['chrom'] == 'chr22'].copy()
        blacklist['same'] = blacklist['start'] // l == blacklist['end'] // l
        blacklist = blacklist[blacklist['same']].copy()
        blacklist['region'] =  blacklist['start'] // l
        blacklist['start'] += blacklist['region'] * s
        blacklist['end'] += blacklist['region'] * s
        blacklist['chrom'] = 'chr23'
        blacklist[['chrom', 'start', 'end']].to_csv(
            output.blacklist, sep='\t', index=False, header=False
        )
        with open(output.mock, 'w') as f:
            # motif shift
            f.write(f'chr23\t{l}\t{l + s}\tecDNA_1_1\t.\t+\n')
            f.write(f'chr23\t{l + s + l}\t{l + s + l + s + 1}\tecDNA_1_2\t.\t-\n')



rule fasta_index:
    output:
        config['workspace'] + '/verification/genome/mock.fa.fai'
    input:
        config['workspace'] + '/verification/genome/mock.fa'
    shell:
        'samtools faidx {input}'


rule bwa_index:
    output:
        config['workspace'] + '/verification/genome/mock.fa.bwt'
    input:
        config['workspace'] + '/verification/genome/mock.fa'
    shell:
        'bwa index {input}'


rule reads_list:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{gsm}_{chrom}.list'
    input:
        bam=rules.merge.output,
        index=rules.index.output
    shell:
        'samtools view -h {input.bam} {wildcards.chrom}'
        ' | samtools sort -O SAM -n -'
        ' | grep -v \'^@\''
        ' | awk \'{{print $1}}\''
        ' | uniq > {output}'


rule extract_reads:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{srr}_{chrom}_{r}.fastq.gz'
    input:
        reads=rules.reads_list.output,
        fq=lambda wildcards: config['samples'][wildcards.gsm][wildcards.srr][f'fq{wildcards.r}']
    params:
        script=lambda wildcards: BASE_DIR + '/tools/extractreads.py',
    shell:
        'python {params.script} {input.fq} {input.reads} {output}'


rule bedgraph:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{gsm}_pileup.bdg'
    input:
        bam=rules.merge.output,
        index=rules.index.output
    shell:
        'bedtools genomecov -bga -ibam {input.bam} > {output}'


rule control_depth:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{gsm}_control_depth.value'
    input:
        rules.bedgraph.output
    shell:
        'grep \'^chr[0-9XY]\{{1,2\}}[[:space:]]\' {input}'
        ' | awk \'{{LENGTH+=$3-$2;BASE+=$4*($3-$2)}}END{{print BASE/LENGTH}}\' > {output}'


rule chrM_depth:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{gsm}_chrM_depth.value'
    input:
        rules.bedgraph.output
    shell:
        'grep \'^chrM\' {input} | awk \'{{LENGTH+=$3-$2;BASE+=$4*($3-$2)}}END{{print BASE/LENGTH}}\' > {output}'


def get_step_num(wildcards):
    filename = config['workspace'] + f'/verification/fastq/{wildcards.prefix}/{wildcards.gsm}/{wildcards.gsm}_control_depth.value'
    with open(filename) as f:
        control_depth = f.read().strip()
    filename = config['workspace'] + f'/verification/fastq/{wildcards.prefix}/{wildcards.gsm}/{wildcards.gsm}_chrM_depth.value'
    with open(filename) as f:
        chrM_depth = f.read().strip()
    ratio = float(control_depth) * float(wildcards.foldchange) / float(chrM_depth)
    return int(1 / ratio) + 1


rule split_reads:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{srr}_{chrom}_{foldchange}_{r}.fastq.gz'
    input:
        fastq=rules.extract_reads.output,
        control=rules.control_depth.output,
        chrM=rules.chrM_depth.output
    params:
        step=get_step_num
    shell:
        'zcat {input.fastq} | awk \'int((NR - 1) / 4) % {params.step} == 0\' | gzip -c > {output}'


rule trim_reads:
    output:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{srr}_{chrom}_L{length}_{r}.fastq.gz'
    input:
        config['workspace'] + '/verification/fastq/{prefix}/{gsm}/{srr}_{chrom}_{r}.fastq.gz'
    shell:
        'zcat {input} | awk \'{{if(NR % 4 % 2 == 0){{print(substr($0, 0, {wildcards.length}))}}else{{print($0)}}}}\' | gzip -c > {output}'


########################################################
rule verification:
    input:
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/calling/{gsm}_ecDNA.bed'
            for gsm in config['samples']
        ],
        [
            config['workspace'] + f'/samples/{gsm[:6]}/{gsm}/verification/{gsm}_ecDNA_deep.value'
            for gsm in config['samples']
        ]


rule extract_tag:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/verification/{gsm}_ecDNA_tag.bed'
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=lambda wildcards: os.path.splitext(config['genome']['fasta'])[0] + '.ecDNA.bed'
    params:
        script=lambda wildcards: BASE_DIR + '/tools/bam2bed.py',
        mapq=config['params']['mapq'],
        include=config['params']['include'],
        exclude=config['params']['exclude'],
        mismatch=config['params']['mismatch']
    shell:
        'python {params.script} -q {params.mapq} -f {params.include} -F {params.exclude} -r {params.mismatch} -L {input.bed} {input.bam} {output}'


rule filter_tag:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/verification/{gsm}_ecDNA_filter.bed'
    input:
        tag=rules.extract_tag.output,
        bed=lambda wildcards: os.path.splitext(config['genome']['fasta'])[0] + '.ecDNA.bed'
    shell:
        'bedtools intersect -a {input.tag} -b {input.bed} > {output}'


rule pileup_tag:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/verification/{gsm}_ecDNA_pileup.bdg'
    input:
        rules.filter_tag.output
    shell:
        'macs2 pileup -i {input} -f BEDPE -o {output}'


rule ecDNA_deep:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/verification/{gsm}_ecDNA_deep.value'
    input:
        rules.pileup_tag.output
    shell:
        'awk \'$4>0{{LENGTH+=$3-$2;BASE+=$4*($3-$2)}}END{{print BASE / LENGTH}}\' {input} > {output}'
