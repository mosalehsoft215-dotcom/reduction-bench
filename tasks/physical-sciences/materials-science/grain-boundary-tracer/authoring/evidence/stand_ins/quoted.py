import json, os, subprocess, sys
answer = sys.argv[2]
subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "real.py"), sys.argv[1], answer])
a = json.load(open(answer))
json.dump({k: str(v) for k, v in a.items()}, open(answer, "w"))
