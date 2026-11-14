#!/bin/bash
# copies the reference relay reduction into place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/setting_recovery.py" /app/solution/setting_recovery.py
