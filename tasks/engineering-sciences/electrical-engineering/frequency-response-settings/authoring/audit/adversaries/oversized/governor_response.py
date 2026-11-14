"""Correct values padded past 32 kB."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref import solve
v = solve(sys.argv[1])
v["pad"] = "x" * 40000
open(sys.argv[2], "w").write(json.dumps(v))
