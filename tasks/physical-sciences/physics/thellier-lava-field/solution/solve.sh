#!/bin/bash
# puts the field reduction in the folder the kept-back lavas are reduced from
set -e
mkdir -p /app/solution
install -m 0644 "$(dirname "$0")/lava_field.py" /app/solution/lava_field.py
