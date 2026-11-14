#!/bin/bash
# copies the reference precharge reduction into place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/charge_audit.py" /app/solution/charge_audit.py
