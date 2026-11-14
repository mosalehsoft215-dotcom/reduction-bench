import os
import sys

target = sys.argv[-1]
os.mkfifo(target)
