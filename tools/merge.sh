#!/bin/bash

set -euxo pipefail

cat $1 $2 $3 | sort | uniq -c | awk '$1=="3" { print $2 }' > $4
