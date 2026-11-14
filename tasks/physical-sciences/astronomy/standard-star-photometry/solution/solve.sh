#!/bin/bash
# Installs the reference reduction.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/photometry_reduce.py" /app/solution/photometry_reduce.py
