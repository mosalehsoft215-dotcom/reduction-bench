#!/bin/bash
# the reference maturity reduction goes to /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/maturity.py" /app/solution/
