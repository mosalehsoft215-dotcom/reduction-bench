#!/bin/bash
# puts the reference uplift reduction in /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/rail_thrust.py" /app/solution/
