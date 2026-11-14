#!/bin/bash
# installs the reference roof tie evaluation
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/tie_loads.py" /app/solution/tie_loads.py
