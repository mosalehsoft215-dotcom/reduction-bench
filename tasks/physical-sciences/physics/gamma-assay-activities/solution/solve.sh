#!/bin/bash
# Install the reference program and run it on the two released assays.
set -e
src=$(dirname "$(readlink -f "$0")")
mkdir -p /app/solution
cp "$src/assay_activity.py" /app/solution/assay_activity.py
for d in GA-04 GA-09; do
  python3 /app/solution/assay_activity.py "/app/data/assays/$d" "/tmp/$d.json" && cat "/tmp/$d.json"
done
