#!/bin/bash
# puts the reference reduction where the sealed tape pieces are run from
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/tape_ic.py" /app/solution/tape_ic.py
