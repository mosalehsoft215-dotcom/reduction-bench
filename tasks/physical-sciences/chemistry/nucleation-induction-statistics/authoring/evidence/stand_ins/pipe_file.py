import json, os, pathlib, sys, time
answer = pathlib.Path(sys.argv[2])
os.mkfifo(answer)
