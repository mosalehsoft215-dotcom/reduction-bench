#!/bin/bash
# installs the reference station reduction
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/notch_flow.py" /app/solution/notch_flow.py
