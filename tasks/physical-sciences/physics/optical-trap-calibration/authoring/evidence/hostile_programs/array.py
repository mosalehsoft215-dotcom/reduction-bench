import pathlib, sys
pathlib.Path(sys.argv[-1]).write_text('[0.05, 0.25, 120.0, 500.0]')
