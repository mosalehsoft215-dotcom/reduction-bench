#!/bin/bash
# installs the reference casting evaluation into /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/cast_maturity.py" /app/solution/
