#!/bin/bash

set -euxo pipefail

awk '{print $4}' $1 | sort | uniq -c | awk '$1=="3" {print $2}' > $2
