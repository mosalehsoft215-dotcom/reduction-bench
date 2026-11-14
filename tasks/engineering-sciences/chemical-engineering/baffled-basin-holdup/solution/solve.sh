#!/bin/bash
# puts the reference basin reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/basin_holdup.py" /app/solution/basin_holdup.py
