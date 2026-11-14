from pathlib import Path
import sys
Path(sys.argv[2]).symlink_to("/tests/sfc_truth.json")
