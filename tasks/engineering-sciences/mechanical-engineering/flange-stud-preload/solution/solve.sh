#!/bin/bash
# installs the reference joint reduction
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/preload_reduce.py" /app/solution/preload_reduce.py
