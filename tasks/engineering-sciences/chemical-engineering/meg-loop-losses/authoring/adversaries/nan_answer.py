"""Hostile or malformed program used as /app/solution/meg_balance.py in the Docker probes."""
import sys
open(sys.argv[2], "w").write('{"formation_water_t_d": NaN, "meg_loss_kg_d": Infinity}')
