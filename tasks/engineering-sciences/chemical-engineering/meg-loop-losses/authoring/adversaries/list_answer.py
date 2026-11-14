"""Hostile or malformed program used as /app/solution/meg_balance.py in the Docker probes."""
import json, sys
json.dump([60.0, 1900.0], open(sys.argv[2], "w"))
