rule accessible_tag:
    output:
        temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_tag.bed')
    input:
        bam=rules.merge.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2bed.py',
        mapq=config['params']['mapq']
    shell:
        'python {params.script} -q {params.mapq} -L {input.bed} {input.bam} {output}'


rule accessible_peak:
    output:
        peak=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks.narrowPeak',
        treat_bdg=temp(config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_treat_pileup.bdg'),
        lambda_bdg=config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_control_lambda.bdg'
    input:
        rules.accessible_tag.output
    params:
        outdir=config['workspace'] + '/samples/{prefix}/{gsm}/accessible',
        name='{gsm}_accessible'
    shell:
        'macs2 callpeak -t {input} -f BEDPE -g hs --keep-dup all --outdir {params.outdir} -n {params.name}'
        ' -B --nomodel -q 0.05 --nolambda --max-gap 200 --min-length 1000'


# def get_accessible_tag_num(wildcards):
#     path = config['workspace'] + f'/samples/{wildcards.prefix}/{wildcards.gsm}/accessible/{wildcards.gsm}_accessible_tag.bed'
#     with open(path, 'rb') as f:
#         count = 0
#         buff = f.read(65536)
#         while buff:
#             count += buff.count('\n')
#             buff = f.read(65536)
#     return count


# rule accessible_rpkm:
#     output:
#         config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_peaks_rpkm.bed'
#     input:
#         peak=rules.accessible_peak.output.peak,
#         chrom_sizes=rules.chrom_sizes.output,
#         tag=rules.accessible_tag.output
#     params:
#         awk=os.path.dirname(workflow.snakefile) + '/tools/rpkm.awk'
#     shell:
#         'bedtools sort -g {input.chrom_sizes} -i {input.peak}'
#         ' | bedtools merge -d 12500 -i stdin -c 4 -o first '
#         ' | bedtools intersect -c -a stdin -b {input.tag}'
#         ' | awk -v TOTAL=`cat {input.tag} | wc -l` -f {params.awk} > {output}'


# rule accessible_super:
#     output:
#         config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_super.bed'
#     input:
#         rules.accessible_rpkm.output
#     params:
#         script=os.path.dirname(workflow.snakefile) + '/tools/superfilter.py'
#     shell:
#         'python {params.script} {input} {output}'


# rule accessible_super:
#     output:
#         config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_super.narrowPeak'
#     input:
#         rules.accessible_peak.output.peak
#     params:
#         script=os.path.dirname(workflow.snakefile) + '/tools/superfilter.py'
#     shell:
#         'python {params.script} {input} {output}'


# rule accessible_merge:
#     output:
#         config['workspace'] + '/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_super.bed'
#     input:
#         peak=rules.accessible_peak.output.peak,
#         chrom_sizes=rules.chrom_sizes.output
#     shell:
#         'bedtools sort -g {input.chrom_sizes} -i {input.peak}'
#         ' | bedtools merge -d 125000 -i stdin -c 4 -o first > {output}'
