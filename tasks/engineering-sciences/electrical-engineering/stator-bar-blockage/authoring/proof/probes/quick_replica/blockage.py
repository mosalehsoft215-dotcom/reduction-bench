import json, sys
import routes
got = routes.main(sys.argv[1], None, routes.ROUTES["quick_replica"])
json.dump({"blocked_bar": got["blocked_bar"], "flow_fraction": float(got["flow_fraction"])}, open(sys.argv[2], "w"))
