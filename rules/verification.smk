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