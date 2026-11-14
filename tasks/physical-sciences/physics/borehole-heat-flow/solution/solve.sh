#!/bin/bash
# Puts the reference borehole reduction where the judge looks for it.
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/well_reduce.py" /app/solution/well_reduce.py
