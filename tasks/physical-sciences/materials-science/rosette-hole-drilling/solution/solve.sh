#!/bin/bash
# Copies the reference reduction into place.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/hole_reduce.py" /app/solution/hole_reduce.py
