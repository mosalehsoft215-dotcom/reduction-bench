#!/bin/bash
set -euo pipefail
mkdir -p /app/solution
cp "$(dirname "$0")/coupon_fit.py" /app/solution/coupon_fit.py
