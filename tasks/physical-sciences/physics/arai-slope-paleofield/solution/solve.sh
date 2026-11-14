#!/bin/bash
# puts the field reduction in the folder the kept-back basalts are reduced from
set -e
mkdir -p /app/solution
install -m 0644 "$(dirname "$0")/arai_field.py" /app/solution/arai_field.py
