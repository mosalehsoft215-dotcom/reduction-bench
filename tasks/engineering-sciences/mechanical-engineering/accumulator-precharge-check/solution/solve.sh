#!/bin/bash
# copies the reference precharge reduction into place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/precharge.py" /app/solution/precharge.py
