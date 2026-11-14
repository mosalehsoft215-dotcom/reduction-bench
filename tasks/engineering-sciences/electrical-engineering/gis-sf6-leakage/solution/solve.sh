#!/bin/bash
# install the bay leak program and try it on the two public bays
set -e
mkdir -p /app/solution /app/output
cp /solution/gis_leak.py /app/solution/
for bay in /app/data/bays/*/; do
  python3 /app/solution/gis_leak.py "$bay" "/app/output/$(basename "$bay").json"
done
