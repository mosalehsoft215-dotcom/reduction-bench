"""Noise redraws: each campaign is filed again with a new noise seed and reduced by both programs.

    python redraws.py 20       -> redraws.txt, redraws.json

The model parameters stay fixed; only the measurement noise changes (seed 33000 + 53 i + campaign index). For every
program, campaign and figure the largest error over the redraws is kept in natural units (kJ/mol for Ea, absolute
for the exponents, a fraction for the delays), so it can be held against the tolerance directly. second_solver.py
needs SciPy.
"""
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAMPAIGNS = ("P1", "P2", "H1", "H2", "H3", "H4", "H5")
FIGURES = ("activation_energy_kj_mol", "pressure_exponent", "oxygen_exponent", "tau_report_1_us", "tau_report_2_us")


def importable():
    for part in ("solution", "authoring/evidence", "authoring/provenance"):
        if str(ROOT / part) not in sys.path:
            sys.path.append(str(ROOT / part))


def miss(answer, truth):
    return [abs(answer[f] / truth[f] - 1) if f.startswith("tau_") else abs(answer[f] - truth[f]) for f in FIGURES]


def one_redraw(campaign, draw):
    importable()
    import file_campaign
    import second_solver
    import shock_reduce
    with tempfile.TemporaryDirectory() as scratch:
        place = Path(scratch)
        truth = file_campaign.file_campaign(campaign, place / "c",
                                            noise_seed=33000 + 53 * draw + CAMPAIGNS.index(campaign))
        first = shock_reduce.reduce(place / "c")
        second_solver.main(str(place / "c"), str(place / "second.json"))
        second = json.loads((place / "second.json").read_text())
    return campaign, miss(first, truth), miss(second, truth)


if __name__ == "__main__":
    draws = int(sys.argv[1])
    largest = {}
    jobs = [(c, d) for c in CAMPAIGNS for d in range(draws)]
    with ProcessPoolExecutor(max_workers=6) as pool:
        for campaign, first, second in pool.map(one_redraw, *zip(*jobs)):
            for program, errs in (("reference", first), ("second", second)):
                for figure, err in zip(FIGURES, errs):
                    slot = (program, campaign, figure)
                    largest[slot] = max(largest.get(slot, 0.0), err)
    report = [f"largest error over {draws} redraws of each campaign "
              "(Ea in kJ/mol, exponents absolute, delays as fractions of the true delay)"]
    for program in ("reference", "second"):
        report += [f"{program} {c}: " + ", ".join(f"{f.split('_')[0]} {largest[(program, c, f)]:.3f}" for f in FIGURES)
                   for c in CAMPAIGNS]
    (HERE / "redraws.txt").write_bytes(("\n".join(report) + "\n").encode())
    flat = {"|".join(slot): value for slot, value in largest.items()}
    (HERE / "redraws.json").write_bytes((json.dumps(flat, indent=1) + "\n").encode())
    print("\n".join(report))
