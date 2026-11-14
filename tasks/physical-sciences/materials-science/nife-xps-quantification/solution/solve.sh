#!/bin/bash
# Puts the reference XPS reduction in place.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/xps_reduce.py" /app/solution/xps_reduce.py
