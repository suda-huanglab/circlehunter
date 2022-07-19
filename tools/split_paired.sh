#!/bin/bash

set -euxo pipefail

grep -wFf $1 $2 | sed 'N;s/\n/\t/' | awk '$1==$10 && $7==$16 && $6>0 && $15>0 { print $4 }' > $3
