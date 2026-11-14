#!/bin/bash
# Installs the reference calibration program where the judge runs it.
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/trap_reduce.py" /app/solution/trap_reduce.py
