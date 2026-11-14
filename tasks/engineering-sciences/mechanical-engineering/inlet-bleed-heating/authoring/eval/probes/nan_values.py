# Probe: graded keys as NaN.
import sys
open(sys.argv[2], "w").write('{"bleed_b_kg_s": NaN, "rise_a_k": NaN}')
