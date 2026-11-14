#!/bin/bash
# copies the reference measurement to where the sealed stations are measured from
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/split_measure.py" /app/solution/split_measure.py
