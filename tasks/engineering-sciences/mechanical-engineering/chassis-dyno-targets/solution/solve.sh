#!/bin/bash
# puts the reference glide reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/dyno_targets.py" /app/solution/dyno_targets.py
