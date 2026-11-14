#!/bin/bash
# installs the reference brine-pulse reduction as the program
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/pulse_sweep.py" /app/solution/pulse_sweep.py
