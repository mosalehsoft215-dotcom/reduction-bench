#!/bin/bash
# Puts the reference reduction where the rig runs it.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/crack_reduce.py" /app/solution/crack_reduce.py
