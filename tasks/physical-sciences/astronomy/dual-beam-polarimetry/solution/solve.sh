#!/bin/bash
# puts the reference reduction where the judge expects the program
set -euo pipefail
install -d /app/solution
install -m 0644 "$(dirname "$0")/pol_reduce.py" /app/solution/pol_reduce.py
