from functools import lru_cache
import os


rule verification_config:
    input:
        config['workspace'] + '/verification/genome/mock.fa',
        config['workspace'] + '/verification/genome/mock.blacklist.bed',
        config['workspace'] + '/verification/genome/mock.fa.fai',
        config['workspace'] + '/verification/genome/mock.fa.bwt',
        [
            config['workspace'] + f'/verification/samples/{gsm[:6]}/{gsm}/fastq/{gsm}_chr22.list'
            for gsm in config['samples']
        ],
        [
            config['workspace'] + f'/verification/samples/{gsm[:6]}/{gsm}/fastq/{srr}_{chrom}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for chrom in ('chr22', 'chrM')
            for r in (1, 2)
        ],
        [
            config['workspace'] + f'/verification/samples/{gsm[:6]}/{gsm}/foldchange/{gsm}_chrM_depth.value'
            for gsm in config['samples']
        ],[
            config['workspace'] + f'/verification/samples/{gsm[:6]}/{gsm}/fastq/{srr}_chrM_{foldchange}_{r}.fastq.gz'
            for gsm in config['samples']
            for srr in config['samples'][gsm]
            for foldchange in range(10, 101, 10)
            for r in (1, 2)
        ]


rule mock_genome:
    output:
        fasta=config['workspace'] + '/verification/genome/mock.fa',
        blacklist=config['workspace'] + '/verification/genome/mock.blacklist.bed'
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
        f.write('>chrF\n')
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
        blacklist[['chrom', 'start', 'end']].to_csv(
            output.blacklist, sep='\t', index=False, header=False
        )



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
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/fastq/{gsm}_{chrom}.list'
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
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/fastq/{srr}_{chrom}_{r}.fastq.gz'
    input:
        reads=rules.reads_list.output,
        fq=lambda wildcards: config['samples'][wildcards.gsm][wildcards.srr][f'fq{wildcards.r}']
    shell:
        'zcat {input.fq}'
        ' | grep -wFf {input.reads} -A 3'
        ' | grep -v \'^--\''
        ' | gzip -c > {output}'


rule chrM_bed:
    output:
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/foldchange/{gsm}_chrM.bed'
    input:
        bam=rules.merge.output,
        index=rules.index.output
    shell:
        'samtools idxstats {input.bam}'
        ' | grep \'chrM\''
        ' | awk \'BEGIN{{OFS="\\t"}}{{print $1, 0, $2}}\' > {output}'


rule chrM_tag:
    output:
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/foldchange/{gsm}_chrM_tag.bed'
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.chrM_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2bed.py',
        mapq=config['params']['mapq'],
        include=config['params']['include'],
        exclude=config['params']['exclude'],
        mismatch=config['params']['mismatch']
    shell:
        'python {params.script} -q {params.mapq} -f {params.include} -F {params.exclude} -r {params.mismatch} -L {input.bed} {input.bam} {output}'


rule chrM_pileup:
    output:
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/foldchange/{gsm}_chrM.bdg'
    input:
        rules.chrM_tag.output
    shell:
        'macs2 pileup -i {input} -f BEDPE -o {output}'


rule chrM_depth:
    output:
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/foldchange/{gsm}_chrM_depth.value'
    input:
        rules.chrM_pileup.output
    shell:
        'awk \'{{LENGTH+=$3-$2;BASE+=$4*($3-$2)}}END{{print BASE / LENGTH}}\' {input} > {output}'


@lru_cache(maxsize=1024)
def line_count(filename, srr):
    srr_re = re.compile(f'^{srr}\.')
    count = 0
    with open(filename) as f:
        for line in f:
            if srr_re.match(line):
                count += 1
    return count


def get_line_num(wildcards):
    import re
    comment_depth = get_accessible_lambda(wildcards)
    filename = config['workspace'] + f'/verification/samples/{wildcards.prefix}/{wildcards.gsm}/foldchange/{wildcards.gsm}_chrM_depth.value'
    with open(filename) as f:
        chrM_depth = f.read().strip()
    ratio = int(wildcards.foldchange) / 100 * float(comment_depth) / float(chrM_depth)
    filename = config['workspace'] + f'/verification/samples/{wildcards.prefix}/{wildcards.gsm}/fastq/{wildcards.gsm}_{wildcards.chrom}.list'
    count = line_count(filename, wildcards.srr)
    return int(count * ratio + 1) * 4


rule split_reads:
    output:
        config['workspace'] + '/verification/samples/{prefix}/{gsm}/fastq/{srr}_{chrom}_{foldchange}_{r}.fastq.gz'
    input:
        rules.extract_reads.output
    params:
        num=get_line_num
    shell:
        'zcat {input} | sed -n \'1,{params.num}p\' | gzip -c > {output}'