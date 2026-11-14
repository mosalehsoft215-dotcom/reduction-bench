"""writes the values as strings and NaN"""
import sys
open(sys.argv[2], "w").write('{"tau0_ns": "3.0", "diffusion_cm2_per_s": NaN, "encounter_radius_A": "6.4"}')
