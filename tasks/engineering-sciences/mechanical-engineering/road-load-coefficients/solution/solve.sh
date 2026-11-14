#!/bin/bash
# installs the reference road-load reduction
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/road_load.py" /app/solution/road_load.py
