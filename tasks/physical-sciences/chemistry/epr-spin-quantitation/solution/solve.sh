#!/bin/bash
# put the spin-count program and its lineshape module where the sessions expect them
mkdir -p /app/solution
cp "$(dirname "$0")/spin_count.py" "$(dirname "$0")/powder.py" /app/solution/
