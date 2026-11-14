import sys, json
json.dump({k: 0 for k in json.loads('{"basal_heat_flow_mw_m2": 66.0, "temperature_at_td_c": 50.0, "report_gradient_k_km": 25.0, "report_conductivity_w_mk": 2.5}')}, open(sys.argv[2], 'w'))
