import pandas as pd

fastq = pd.read_csv('/home/andy/Downloads/ATACdb/fq.txt', names=['gsm', 'srr', 'ftp', 'md5'], comment='#')
fastq['r'] = fastq['ftp'].str.extract(r'(\d).fastq.gz')
fastq = fastq.set_index(['srr', 'r'])

root = '/mnt/2w/data1/andy/raw_data/ncbi'


rule download_all:
    input:
        [
            f'{root}/{row["gsm"][:6]}/{row["gsm"]}/{srr}/{srr}_{r}.fastq.gz'
            for (srr, r), row in fastq.iterrows()
        ]


rule download_fastq:
    output:
        root + '/{prefix}/{gsm}/{srr}/{srr}_{r}.fastq.gz'
    params:
        ftp=lambda wildcards: 'ftp://' + fastq.loc[(wildcards.srr, wildcards.r), 'ftp'],
        md5=lambda wildcards: fastq.loc[(wildcards.srr, wildcards.r), 'md5'],
        tmp=root + '/{prefix}/{gsm}/{srr}/{srr}_{r}.fastq.gz.tmp',
        checksum=root + '/{prefix}/{gsm}/{srr}/{srr}_{r}.fastq.gz.md5'
    shell:
        'wget -c -w 30 -T 60 -O {params.tmp} {params.ftp}'
        ' && echo "{params.md5}  {params.tmp}" > {params.checksum}'
        ' && md5sum -c {params.checksum}'
        ' && mv {params.tmp} {output}'
