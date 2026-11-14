"""Regenerates the seven runs from their seeds into a scratch folder and compares every file with the bundle's
copy by sha256; writes regenerate_check.txt.

    python regenerate_check.py
"""
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
GEN = BUNDLE / "authoring" / "calorimeter" / "campaign.py"


def bundle_folder(run):
    p = BUNDLE / "environment" / "data" / "runs" / f"run_{run}"
    return p if p.exists() else BUNDLE / "tests" / "runs" / run


lines = []
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run([sys.executable, str(GEN), "--redraw", "0", "--out", tmp], check=True, cwd=GEN.parent)
    for fresh in sorted(Path(tmp).iterdir()):
        for f in sorted(fresh.iterdir()):
            kept = bundle_folder(fresh.name) / f.name
            a = hashlib.sha256(f.read_bytes()).hexdigest()
            b = hashlib.sha256(kept.read_bytes()).hexdigest() if kept.exists() else "missing"
            lines.append(f"{'same' if a == b else 'DIFFERENT'}  {a[:16]}  {kept.relative_to(BUNDLE).as_posix()}")
same = sum(x.startswith("same") for x in lines)
lines.append(f"{same} of {len(lines)} files reproduced byte for byte")
(HERE / "regenerate_check.txt").write_text("\n".join(lines) + "\n", newline="\n")
print("\n".join(lines))
