#!/bin/bash
# Installs the reference reduction.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/resonator_reduce.py" /app/solution/resonator_reduce.py
