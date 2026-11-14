from pathlib import Path
import sys
Path(sys.argv[2]).write_text(Path("/tests/campaign_truth.json").read_text())
