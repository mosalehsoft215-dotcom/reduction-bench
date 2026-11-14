#!/bin/bash
# installs the reference lift reduction in /app/solution
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/stress_free.py" /app/solution/
