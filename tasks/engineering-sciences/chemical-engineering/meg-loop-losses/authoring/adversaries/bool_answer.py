"""Hostile or malformed program used as /app/solution/meg_balance.py in the Docker probes."""
import json, sys
json.dump({"formation_water_t_d": True, "meg_loss_kg_d": True}, open(sys.argv[2], "w"))
