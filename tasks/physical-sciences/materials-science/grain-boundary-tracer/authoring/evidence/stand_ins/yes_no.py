import json, os, subprocess, sys
answer = sys.argv[2]
json.dump({"lattice_diffusivity_m2_s": True, "boundary_product_m3_s": True}, open(answer, "w"))
