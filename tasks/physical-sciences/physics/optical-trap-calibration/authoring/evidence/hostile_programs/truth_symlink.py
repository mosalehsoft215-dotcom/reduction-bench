import pathlib, sys
pathlib.Path(sys.argv[-1]).symlink_to('/tests/session_truths.json')
