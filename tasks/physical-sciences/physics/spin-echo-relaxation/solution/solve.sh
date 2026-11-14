#!/bin/bash
# puts the reference reduction where the withheld measurements are reduced from
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/echo_reduce.py" /app/solution/echo_reduce.py
