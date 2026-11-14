#!/bin/bash
# copies the reference precharge reduction into place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/accu_gas.py" /app/solution/accu_gas.py
