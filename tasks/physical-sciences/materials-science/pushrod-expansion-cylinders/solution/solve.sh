#!/bin/bash
# Install the CTE reduction and show its results on the two worked cycles.
set -e
d="$(cd "$(dirname "$0")" && pwd)"
mkdir -p /app/solution
cp "$d/cte_reduce.py" /app/solution/cte_reduce.py
for r in DL-06 DL-10; do
  python3 /app/solution/cte_reduce.py "/app/data/cycles/$r" "/tmp/$r.json"
  cat "/tmp/$r.json"
done
