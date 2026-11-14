#!/bin/bash
# Installs the reference reduction.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/ew_reduce.py" /app/solution/ew_reduce.py
