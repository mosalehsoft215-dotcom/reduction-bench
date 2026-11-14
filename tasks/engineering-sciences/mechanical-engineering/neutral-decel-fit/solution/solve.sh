#!/bin/bash
# puts the reference deceleration reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/decel_fit.py" /app/solution/decel_fit.py
