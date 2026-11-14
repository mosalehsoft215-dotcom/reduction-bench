import os, sys
os.mkfifo(sys.argv[2], 0o644)
