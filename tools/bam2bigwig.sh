#!/usr/bin/env bash
set -euxo pipefail

scale=$(samtools idxstats $1 | awk '{sum+=$3}END{print sum/1000000}')

samtools idxstats $1 | awk 'BEGIN{OFS="\t"}{print $1, $2}' > "$2.chrom_sizes"

bedtools genomecov -ibam $1 -bg -scale $scale > "$2.bdg"

bedSort "$2.bdg" "$2.sorted.bdg"

bedGraphToBigWig "$2.sorted.bdg" "$2.chrom_sizes" $2

rm "$2.chrom_sizes" "$2.bdg" "$2.sorted.bdg"
