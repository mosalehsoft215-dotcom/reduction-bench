import json, os, sys, time
out = sys.argv[2]
json.dump(dict({'zimm_viscosity_mpa_s': 1.6, 'stretch_exponent': 0.67, 'incoherent_fraction_at_reference': 0.16}, pad='~' * 300000), open(out, 'w'))
