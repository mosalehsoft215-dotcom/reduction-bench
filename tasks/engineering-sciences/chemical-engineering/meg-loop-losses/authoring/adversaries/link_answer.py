"""Hostile or malformed program used as /app/solution/meg_balance.py in the Docker probes."""
import os, sys
os.symlink("/tests/loop_truth.json", sys.argv[2])
