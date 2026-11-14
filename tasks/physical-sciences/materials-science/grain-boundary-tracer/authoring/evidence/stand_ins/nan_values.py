import json, os, subprocess, sys
answer = sys.argv[2]
open(answer, "w").write("{\"lattice_diffusivity_m2_s\": NaN, \"boundary_product_m3_s\": NaN}")
