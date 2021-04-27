rule cleanup:
    output:
        config['workspace'] + '/samples/{prefix}/{gsm}/cleanup.list'
    input:
        rule=rules.annotate.output
    run:
        with open(output[0], 'w') as f:
            prefix = wildcards.gsm[:6]
            gsm = wildcards.gsm
            for srr in config['samples'][gsm]:
                # trimed fastq
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/preprocess/{srr}_r1.trimed.fq.gz\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/preprocess/{srr}_r2.trimed.fq.gz\n')
                # srr bam
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/mapping/{srr}.sorted.bam\n')
                # clean bam
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/mapping/{gsm}_clean.sorted.bam\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/mapping/{gsm}_clean.sorted.bam.bai\n')
                # accessible
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_tag.bed\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_treat_pileup.bdg\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/accessible/{gsm}_accessible_control_lambda.bdg\n')
                # largeinsert
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_tag.bed\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_pileup.bdg\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_ratio.bdg\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_lambda.bdg\n')
                f.write(config['workspace'] + f'/samples/{prefix}/{gsm}/largeinsert/{gsm}_largeinsert_qvalue.bdg\n')
