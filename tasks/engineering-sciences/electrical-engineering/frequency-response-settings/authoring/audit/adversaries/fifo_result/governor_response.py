"""Leaves a named pipe where the result file should be."""
import os, sys
os.mkfifo(sys.argv[2])
