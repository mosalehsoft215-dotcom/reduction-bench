"""Exit 0 only if the CTRF summary and the JUnit report both show exactly N tests, all passed.

    python3 count_reports.py <ctrf.json> <junit.xml> N

CTRF merges parametrised cases under one name, so the per-case count is taken from JUnit as well.
"""
import json
import sys
from xml.etree import ElementTree


def ctrf_ok(path, n):
    with open(path) as fh:
        summary = json.load(fh)["results"]["summary"]
    quiet = all(summary.get(k, 0) == 0 for k in ("failed", "skipped", "pending", "other"))
    return quiet and summary.get("passed") == summary.get("tests") and summary.get("tests", 0) >= 1


def junit_ok(path, n):
    root = ElementTree.parse(path).getroot()
    suites = [root] if root.tag == "testsuite" else root.findall("testsuite")
    total = sum(int(s.get("tests", 0)) for s in suites)
    bad = sum(int(s.get(k, 0)) for s in suites for k in ("failures", "errors", "skipped"))
    return total == n and bad == 0


if __name__ == "__main__":
    want = int(sys.argv[3])
    ok_c, ok_j = ctrf_ok(sys.argv[1], want), junit_ok(sys.argv[2], want)
    print(f"ctrf {'agrees' if ok_c else 'disagrees'}; junit {'agrees' if ok_j else 'disagrees'} (want {want})")
    sys.exit(0 if ok_c and ok_j else 1)
