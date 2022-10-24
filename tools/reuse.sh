#!/usr/bin/env bash
set -euxo pipefail
# reuse all ecDNA simulation files generate from last run
# when change ecDNA mock number but not seed, you don't need to regenerate all files
# touch files in this sequential order to fool snakemake

find simulation/ecDNA/mock_ecDNA.fa -name 'ecDNA_*.fa' | xargs -r touch
find simulation/ecDNA_fastq -name '*.fq.gz' | xargs -r touch
if [ -d "simulation/ecDNA_preprocess" ]
then
    find simulation/ecDNA_preprocess -name '*.fq.gz' | xargs -r touch
fi
find simulation/ecDNA_qc -name '*.json' | xargs -r touch
find simulation/ecDNA_qc -name '*.html' | xargs -r touch
find simulation/ecDNA_mapping -name '*.sorted.bam' | xargs -r touch
