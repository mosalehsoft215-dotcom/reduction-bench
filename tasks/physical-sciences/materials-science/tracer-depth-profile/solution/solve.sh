#!/bin/bash
# copies the evaluation script to the place later coupons are evaluated from
set -e
mkdir -p /app/solution
cp "$(dirname "$0")/depth_profile.py" /app/solution/depth_profile.py
