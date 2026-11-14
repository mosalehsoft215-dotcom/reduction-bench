#!/bin/bash
# puts the reference free-roll reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/freeroll_fit.py" /app/solution/freeroll_fit.py
