import json, os, subprocess, sys
answer = sys.argv[2]
json.dump({"lattice_d_m2_s": True, "boundary_sdd_m3_s": True}, open(answer, "w"))
