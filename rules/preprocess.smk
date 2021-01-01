rule trim:
    output:
        fq1=config['workspace'] + '/samples/{prefix}/{gsm}/preprocess/{srr}_r1.trimed.fq.gz',
        fq2=config['workspace'] + '/samples/{prefix}/{gsm}/preprocess/{srr}_r2.trimed.fq.gz',
        json=config['workspace'] + '/samples/{prefix}/{gsm}/qc/{srr}_fastp.json',
        html=config['workspace'] + '/samples/{prefix}/{gsm}/qc/{srr}_fastp.html'
    input:
        adapter=config['adapter'],
        fq1=lambda wildcards: config['samples'][wildcards.gsm][wildcards.srr]['fq1'],
        fq2=lambda wildcards: config['samples'][wildcards.gsm][wildcards.srr]['fq2']
    log:
        config['workspace'] + '/samples/{prefix}/{gsm}/log/{srr}_fastp.log'
    threads: 8 if workflow.cores > 8 else workflow.cores
    shell:
        'fastp -g -q 5 -u 50 -n 15 -w {threads} --adapter_fasta {input.adapter}'
        ' -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}'
        ' -j {output.json} -h {output.html} 2>{log}'
