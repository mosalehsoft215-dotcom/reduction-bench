#!/bin/bash
# Puts the reference reduction where the held-back anneals will look for it.
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/profile_fit.py" /app/solution/profile_fit.py
