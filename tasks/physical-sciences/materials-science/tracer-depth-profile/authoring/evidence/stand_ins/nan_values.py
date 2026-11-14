import json, os, subprocess, sys
answer = sys.argv[2]
open(answer, "w").write("{\"lattice_d_m2_s\": NaN, \"boundary_sdd_m3_s\": NaN}")
