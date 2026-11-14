import json, os, sys, time
out = sys.argv[2]
if os.fork() == 0:
    os.setsid()
    time.sleep(3000)
    os._exit(0)
