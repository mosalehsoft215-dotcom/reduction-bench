"""Feeds both judges (test_series.py through pytest, verdict.py) a table of malformed and near-miss result sets.
Run as root in the verifier image:  python3 broken_answers.py > broken_answers.csv
Only the row "true values" should be accepted, and only when both judges accept it.
"""
import json
import math
import shutil
import subprocess
from pathlib import Path

PICKED = Path("/var/lib/eff-check/results")
SPEC = json.loads(Path("/tests/effusion_truth.json").read_text())
H, P = "dH_sub_298_kJ_mol", "p_298_Pa"
LH, LP = SPEC["limits_percent"][H], SPEC["limits_percent"][P]


def obj(h, p):
    return json.dumps({H: h, P: p})


# (label, function of the true pair giving the file text, series to leave out)
TABLE = [
    ("true values", lambda h, p: obj(h, p), None),
    ("dH 1.03 limits off", lambda h, p: obj(h * (1 + 1.03 * LH / 100), p), None),
    ("p 1.03 limits off", lambda h, p: obj(h, p * (1 - 1.03 * LP / 100)), None),
    ("dH in J/mol", lambda h, p: obj(h * 1000, p), None),
    ("ln p instead of p", lambda h, p: obj(h, math.log(p)), None),
    ("p in kPa", lambda h, p: obj(h, p / 1000), None),
    ("one series without result", lambda h, p: obj(h, p), "h1"),
    ("NaN", lambda h, p: '{"%s": NaN, "%s": NaN}' % (H, P), None),
    ("Infinity", lambda h, p: '{"%s": Infinity, "%s": 0.1}' % (H, P), None),
    ("numbers as text", lambda h, p: json.dumps({H: str(h), P: str(p)}), None),
    ("true/false", lambda h, p: '{"%s": true, "%s": false}' % (H, P), None),
    ("key twice", lambda h, p: '{"%s": 1, "%s": %r, "%s": %r}' % (H, H, h, P, p), None),
    ("list", lambda h, p: json.dumps([h, p]), None),
    ("empty file", lambda h, p: "", None),
    ("zeros", lambda h, p: obj(0, 0), None),
    ("signs flipped", lambda h, p: obj(-h, -p), None),
]

print("case,pytest_exit,verdict_exit,refused_by_both")
for label, make, skip in TABLE:
    shutil.rmtree(PICKED, ignore_errors=True)
    PICKED.mkdir(parents=True)
    for series, true in SPEC["series"].items():
        if series != skip:
            (PICKED / f"{series}.json").write_text(make(true[H], true[P]))
    first = subprocess.run(["python3", "-B", "-I", "-m", "pytest", "-q", "--noconftest", "-p", "no:cacheprovider",
                            "/tests/test_series.py"], capture_output=True).returncode
    second = subprocess.run(["python3", "-B", "-I", "/tests/verdict.py"], capture_output=True).returncode
    print(f"{label},{first},{second},{'yes' if first and second else 'no'}")
