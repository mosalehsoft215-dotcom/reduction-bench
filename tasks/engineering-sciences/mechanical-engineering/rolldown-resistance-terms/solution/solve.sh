#!/bin/bash
# puts the reference rolldown reduction in place
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/rolldown_fit.py" /app/solution/rolldown_fit.py
