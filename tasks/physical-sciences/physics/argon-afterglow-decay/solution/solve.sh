#!/bin/bash
set -eu
folder="$(cd "$(dirname "$0")" && pwd)"
install -d /app/solution
install -m 644 "$folder/decay_fit.py" /app/solution/decay_fit.py
