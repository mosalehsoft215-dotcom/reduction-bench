#!/bin/bash
# installs the reference member force evaluation
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/canopy_tension.py" /app/solution/canopy_tension.py
cp "$(dirname "$0")/held_table.npy" /app/solution/held_table.npy
