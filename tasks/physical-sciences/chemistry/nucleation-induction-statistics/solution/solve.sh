#!/bin/bash
# installs the reference reduction where the held-back runs will look for it
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/induction_fit.py" /app/solution/induction_fit.py
