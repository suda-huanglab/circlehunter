import os


rule all_o:
    input:
        expand(
            config['workspace'] + '/samples/{sample}/calling/{sample}.microDNA.bed',
            sample=config['samples']
        )


rule trim_o:
    output:
        fq1=config['workspace'] + '/samples/{sample}/preprocess/{sample}_r1.trimed.fq.gz',
        fq2=config['workspace'] + '/samples/{sample}/preprocess/{sample}_r2.trimed.fq.gz',
        json=config['workspace'] + '/samples/{sample}/qc/{sample}_fastp.json',
        html=config['workspace'] + '/samples/{sample}/qc/{sample}_fastp.html'
    input:
        adapter=config['adapter'],
        fq1=lambda wildcards: config['samples'][wildcards.sample]['fq1'],
        fq2=lambda wildcards: config['samples'][wildcards.sample]['fq2']
    log:
        config['workspace'] + '/samples/{sample}/log/fastp.log'
    threads: 8
    shell:
        'fastp -g -q 5 -u 50 -n 15 -w {threads} --adapter_fasta {input.adapter}'
        ' -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}'
        ' -j {output.json} -h {output.html} 2>{log}'


rule mapping_o:
    output:
        map=pipe(config['workspace'] + '/samples/{sample}/mapping/{sample}.map.sam'),
        split=pipe(config['workspace'] + '/samples/{sample}/mapping/{sample}.split.sam')
    input:
        fq1=rules.trim_o.output.fq1,
        fq2=rules.trim_o.output.fq2
    log:
       bwa=config['workspace'] + '/samples/{sample}/log/bwa.log',
       samblaster=config['workspace'] + '/samples/{sample}/log/samblaster.log'
    params:
        rg='\'@RG\\tID:{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:ILLUMINA\'',
        index=config['genomes']['hg38']['bwa_index']
    threads: 30
    priority: -1
    shell:
        'bwa mem -t {threads} -R {params.rg} {params.index} {input.fq1} {input.fq2} 2>{log.bwa}'
        ' | samblaster -e --minNonOverlap 10 -s {output.split} > {output.map} 2>{log.samblaster}'


rule sam2bam_o:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.{category}.bam'
    input:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.{category}.sam'
    shell:
        'samtools view -bS {input} > {output}'


rule clean_bam_o:
    output:
          config['workspace'] + '/samples/{sample}/mapping/{sample}.clean.bam'
    input:
         config['workspace'] + '/samples/{sample}/mapping/{sample}.map.bam'
    threads: 8
    shell:
         'samtools view -@ {threads} -b -f 2 {input} > {output}'


rule split_bed_o:
    output:
        config['workspace'] + '/samples/{sample}/calling/{sample}.split.bed'
    input:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.split.bam'
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/split_bed.sh'
    shell:
         'bash {params.script} {input} {output}'


rule concordant_bed_o:
    output:
          config['workspace'] + '/samples/{sample}/calling/{sample}.concordant.bed'
    input:
        rules.clean_bam_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/concordant_bed.sh'
    shell:
         'bash {params.script} {input} {output}'


rule split_freq2_o:
    output:
        config['workspace'] + '/samples/{sample}/calling/{sample}.split.list'
    input:
        rules.split_bed_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/split_freq2.sh',
    shell:
        'bash {params.script} {input} {output}'


rule concordant_freq3_o:
    output:
        config['workspace'] + '/samples/{sample}/calling/{sample}.concordant.list'
    input:
        rules.concordant_bed_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/concordant_freq3.sh',
    shell:
        'bash {params.script} {input} {output}'


rule split_paired_o:
    output:
        config['workspace'] + '/samples/{sample}/calling/{sample}.paired.list'
    input:
        rules.split_freq2_o.output,
        rules.split_bed_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/split_paired.sh',
    shell:
        'bash {params.script} {input} {output}'


rule merge_o:
    output:
        config['workspace'] + '/samples/{sample}/calling/{sample}.merge.list'
    input:
        rules.split_freq2_o.output,
        rules.concordant_freq3_o.output,
        rules.split_paired_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/merge.sh',
    shell:
         'bash {params.script} {input} {output}'


rule calling_o:
    output:
          config['workspace'] + '/samples/{sample}/calling/{sample}.microDNA.bed'
    input:
         rules.merge_o.output,
         rules.concordant_bed_o.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/calling.sh',
    shell:
         'bash {params.script} {input} {output}'

