#!/bin/bash
# Installs the reference vault reduction.
set -e
mkdir -p /app/solution
cp "$(dirname "$(readlink -f "$0")")/vault_reduce.py" /app/solution/vault_reduce.py
