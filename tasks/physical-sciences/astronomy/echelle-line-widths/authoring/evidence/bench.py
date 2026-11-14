"""Error budget of the echelle task in bounds (0.3 km/s on the velocity, 2 % on the width; per exposure the larger)
for the reference, the second reduction and the shortcuts, plus the worst case over eight noise redraws (six worker
processes).  Writes bench.json next to this file.  Usage: python bench.py"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "solution"), str(HERE)]
import ew_reduce as ref  # noqa: E402
from second_ew import reduce_exposure as second  # noqa: E402

PUBLIC = ("S-04", "S-09")


def layout(base):
    truth = json.loads((base / "environment/data/spectra/worked_values.json").read_text())
    truth.update(json.loads((base / "tests/spectrum_truth.json").read_text()))
    return truth, {n: base / ("environment/data/spectra" if n in PUBLIC else "tests/sealed_spectra") / n for n in truth}


def miss(got, want):
    return max(abs(got["radial_velocity_kms"] - want["radial_velocity_kms"]) / 0.3,
               abs(got["ew_7665_ma"] / want["ew_7665_ma"] - 1) / 0.02)


def divided(folder, scale):
    """The recorded standard divided out (raised to the airmass ratio when scale), with the offsets, scattered light
    and continua of the full fits, then the stellar lines fitted on the corrected, normalised target."""
    head = json.loads((folder / "exposure.json").read_text())
    rest = json.loads((folder / "lines.json").read_text())["stellar_rest_air_a"]
    ps, pt = ref.reduce_exposure(folder, full=True)
    lam_s, cnt_s = ref.read_order(folder / "standard_order.csv")
    lam_t, cnt_t = ref.read_order(folder / "target_order.csv")
    norm_s = (cnt_s - ps[1]) / ref.bin_pix(ref.cubic(ref.fine_grid(lam_s, ps[0]), ps[2:6]))
    t_s = np.interp(lam_t - pt[0], lam_s - ps[0], norm_s)
    power = head["target_airmass"] / head["standard_airmass"] if scale else 1.0
    t_s = np.clip(t_s, 0.03, None) ** power
    flat = (cnt_t - pt[1]) / ref.bin_pix(ref.cubic(ref.fine_grid(lam_t, pt[0]), pt[2:6])) / t_s
    true_l = lam_t - pt[0]
    use = (true_l > 7635) & (t_s > 0.08)

    def resid(q):
        rv, b, *tau = q
        s = 1 + rv / ref.C_KMS
        model = np.exp(-sum(t * np.exp(-((true_l - c * s) / b) ** 2) for c, t in zip(rest, tau)))
        return (model - flat)[use]
    sol = least_squares(resid, np.r_[pt[8], pt[9], np.full(len(rest), 0.4)],
                        bounds=(np.r_[-150, 0.03, np.zeros(len(rest))], np.r_[150, 0.4, np.full(len(rest), 20)]))
    rv = sol.x[0]
    centre = rest[0] * (1 + rv / ref.C_KMS)
    win = np.abs(true_l - centre) <= 0.9
    return {"radial_velocity_kms": float(rv), "ew_7665_ma": float(np.sum(1 - flat[win]) * ref.PIX * 1e3)}


ROUTES = {
    "the laws as the worked exposures show them: one air temperature for both exposures, no gas toward the target":
        lambda f: ref.reduce_exposure(f, temp=False, ism_t=False),
    "the standard's O2 depths carried to the target unchanged (one air temperature)":
        lambda f: ref.reduce_exposure(f, temp=False),
    "target taken as free of interstellar K I": lambda f: ref.reduce_exposure(f, ism_t=False),
    "one interstellar cloud toward the target (a single Gaussian doublet)": lambda f: ref.reduce_exposure(f, clouds=1),
    "every other step right, no moonlight, fringes or gas toward the target modelled":
        lambda f: ref.reduce_exposure(f, extras=False),
    "the notes' laws with the header airmass, R 60 000, one air temperature and featureless stars":
        lambda f: ref.reduce_exposure(f, free_power=False, free_close=False, ism=False, temp=False, ism_t=False),
    "standard taken as featureless (no interstellar K I fitted)": lambda f: ref.reduce_exposure(f, ism=False),
    "target airmass held at the header value": lambda f: ref.reduce_exposure(f, free_close=False),
    "one constant fitted airmass for the target (no time mean)": lambda f: ref.reduce_exposure(f, ramp=False),
    "resolving power fixed at the evenly filled slit's 60 000": lambda f: ref.reduce_exposure(f, free_power=False),
    "standard's resolving power used for the target": lambda f: ref.reduce_exposure(f, shared_power=True),
    "recorded standard divided out as it is": lambda f: divided(f, False),
    "recorded standard raised to the airmass ratio, then divided out": lambda f: divided(f, True),
}


def job(task):
    kind, folder = task[0], Path(task[1])
    if kind == "reference":
        return ref.reduce_exposure(folder)
    if kind == "second":
        return second(folder)
    return ROUTES[kind](folder)


def main():
    from concurrent.futures import ProcessPoolExecutor
    truth, where = layout(ROOT)
    names = sorted(truth, key=lambda n: (n not in PUBLIC, n))
    with tempfile.TemporaryDirectory() as tmp:
        draws = []
        for seed in range(1, 9):
            out = Path(tmp) / f"draw{seed}"
            subprocess.run([sys.executable, str(ROOT / "authoring/provenance/make_spectra.py"), "--seed", str(seed),
                            "--out", str(out)], check=True)
            draws.append(layout(out))
        tasks = [(k, str(where[n])) for k in ("reference", "second", *ROUTES) for n in names]
        tasks += [(k, str(w2[n])) for t2, w2 in draws for k in ("reference", "second") for n in t2]
        with ProcessPoolExecutor(max_workers=6) as pool:
            got = dict(zip(tasks, pool.map(job, tasks)))
    doc = {"shipped": {n: {k: round(miss(got[(k, str(where[n]))], truth[n]), 3) for k in ("reference", "second")}
                       for n in names},
           "routes": {lab: {n: round(miss(got[(lab, str(where[n]))], truth[n]), 2) for n in names} for lab in ROUTES}}
    sealed = [n for n in names if n not in PUBLIC]
    doc["redraws"] = {k: round(max(miss(got[(k, str(w2[n]))], t2[n]) for t2, w2 in draws for n in t2), 3)
                      for k in ("reference", "second")}
    doc["worked_copy_passes"] = int(sum(miss(truth[p], truth[s]) <= 1 for p in PUBLIC for s in sealed))
    (HERE / "bench.json").write_text(json.dumps(doc, indent=1) + "\n", newline="\n")
    for lab, row in doc["routes"].items():
        print(lab[:62].ljust(62), " ".join(f"{row[n]:6.2f}" for n in names))
    print("shipped", doc["shipped"])
    print("redraws", doc["redraws"], "worked copies", doc["worked_copy_passes"])
    print("per exposure worst redraw", {k: {n: round(max(miss(got[(k, str(w2[n]))], t2[n]) for t2, w2 in draws), 2)
                                            for n in names} for k in ("reference", "second")})


if __name__ == "__main__":
    main()
