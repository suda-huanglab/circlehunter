#!/bin/bash

set -euxo pipefail

bedtools bamtobed -cigar -i $1 \
| sed -e s/_2\\/2/\ 2/g \
| sed -e s/_1\\/1/\ 1/g \
| awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' \
| awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", $8)} 1' \
| awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", $8)} 1' \
| awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", $8)} 1' \
| awk 'BEGIN{FS=OFS=" "} {if (($9=="M" && $NF=="H") || ($9=="M" && $NF=="S"))  {printf ("%s\tfirst\n",$0)} else if (($9=="S" && $NF=="M") || ($9=="H" && $NF=="M")) {printf ("%s\tsecond\n",$0)} }' \
| awk 'BEGIN{FS=OFS="\t"} {gsub("\ ", "", $8)} 1' > $2
