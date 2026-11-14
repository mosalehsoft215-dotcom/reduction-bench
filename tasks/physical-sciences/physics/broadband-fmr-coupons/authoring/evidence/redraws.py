"""Noise redraws: the generator is rerun with other seeds (same films, same sessions) into a scratch folder and
both reductions are scored on every coupon, in allowances (alpha 3 % relative, g 0.0008, mu0 Meff 1.0 mT).

    python redraws.py <scratch dir> [seed offsets...]      -> redraws.txt beside this script
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "solution"))
sys.path.insert(0, str(HERE))
import fmr_reduce as R  # noqa: E402
from routes import errors  # noqa: E402
from second_reduction import film_constants  # noqa: E402


def main():
    scratch = Path(sys.argv[1])
    seeds = [int(s) for s in sys.argv[2:]] or [1, 2, 3]
    lines = []
    for seed in seeds:
        out = scratch / f"seed{seed}"
        subprocess.run([sys.executable, str(ROOT / "authoring/provenance/make_coupons.py"), str(seed), str(out)],
                       check=True)
        truth = json.loads((out / "truth.json").read_text())
        for name, tr in truth.items():
            for label, fn in (("reference", lambda f: R.reduce(f)[0]), ("second", film_constants)):
                e = errors(fn(str(out / name)), tr)
                k = max(e, key=e.get)
                lines.append(f"seed {seed} {name} {label:9s} worst {e[k]:.2f} ({k})")
                print(lines[-1], flush=True)
    (HERE / "redraws.txt").write_text("\n".join(lines) + "\n", newline="\n")


if __name__ == "__main__":
    main()
