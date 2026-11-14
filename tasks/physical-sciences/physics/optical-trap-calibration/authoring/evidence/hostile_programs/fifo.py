import os, sys
os.mkfifo(sys.argv[-1], 0o600)
