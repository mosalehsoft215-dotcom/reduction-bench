#!/bin/bash
# puts the reference jacking reduction in /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/rail_load.py" /app/solution/
