#!/bin/bash
# Installs the reference reduction where the track desk expects it.
mkdir -p /app/solution
cp "$(dirname "$0")/track_reduce.py" /app/solution/track_reduce.py
