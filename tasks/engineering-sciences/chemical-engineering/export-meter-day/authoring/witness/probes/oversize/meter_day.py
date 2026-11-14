import sys
with open(sys.argv[2], "w") as fh:
    fh.write('{"volume_run_a_ksm3": 8000.0, "pad": "' + "x" * 200000 + '"}')
