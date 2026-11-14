"""writes each key twice"""
import sys
open(sys.argv[2], "w").write('{"tau0_ns": 1, "tau0_ns": 3.0, "diffusion_cm2_per_s": 1e-5, "diffusion_cm2_per_s": 1.2e-5,'
                             ' "encounter_radius_A": 1, "encounter_radius_A": 6.4}')
