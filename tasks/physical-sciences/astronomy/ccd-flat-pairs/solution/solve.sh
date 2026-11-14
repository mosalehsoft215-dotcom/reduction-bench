#!/bin/bash
# Installs the reference reduction where the sealed sensors are reduced from.
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/bf_reduce.py" /app/solution/bf_reduce.py
