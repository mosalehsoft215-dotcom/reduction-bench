#!/bin/bash
# puts the reference valve reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/kv_fl_reduce.py" /app/solution/kv_fl_reduce.py
