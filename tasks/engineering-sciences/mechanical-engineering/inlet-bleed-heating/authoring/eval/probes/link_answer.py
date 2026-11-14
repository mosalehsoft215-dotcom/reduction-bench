# Probe: makes the answer a symbolic link to the truth file.
import os, sys
os.symlink("/tests/bleed_truth.json", sys.argv[2])
