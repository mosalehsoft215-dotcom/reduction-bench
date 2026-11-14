#!/bin/bash
# installs the reference casting evaluation into /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/strip_time.py" /app/solution/
