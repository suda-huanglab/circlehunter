#!/bin/bash

set -euxo pipefail

bedtools bamtobed -cigar -i $1 \
| sed -e s/\\//\ /g \
| awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > $2
