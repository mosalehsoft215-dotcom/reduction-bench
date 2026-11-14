"""Readings that pass the two examples, measured on every loop against the truth (in limits).

    python misreadings.py <loops root ...>    # folders holding loop-xx sub-folders; truths from the json files

Each reading changes one decision of the reference and nothing else.
"""
import copy
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "solution"))
import meg_balance as R  # noqa: E402

LIM = {"formation_water_t_d": 0.03, "meg_loss_kg_d": 0.05}


def inventory_route(L, purity=0.995):
    """MEG loss as make-up less the change of the MEG held in the tanks and the two lines."""
    sh, x = L.sheet, L.x
    dst = sh["formation_water_analysis"]
    cl_l, ca_l = L.lab_series("lean MEG", "chloride", 1e-6), L.lab_series("lean MEG", "calcium", 1e-6)
    tds_l = R.tds_from(cl_l, ca_l, dst)
    wm_l = R.w_meg_from_rho(x["DT-101"], x["TT-101"], tds_l)
    cl_r, ca_r = L.lab_series("rich MEG", "chloride", 1e-6), L.lab_series("rich MEG", "calcium", 1e-6)
    tds_r = R.tds_from(cl_r, ca_r, dst)
    wm_r = R.w_meg_from_rho(x["DT-201"], x["TT-201"], tds_r)
    wm_tr = R.w_meg_from_rho(x["DT-401"], x["TT-202"], tds_r)
    lean, rich = x["FI-101"], x["FI-201"] * x["DT-201"]
    vs = L.cum(lean / R.rho(R.T_LINE, wm_l, tds_l))
    vf = L.cum(rich / R.rho(R.T_LINE, wm_r, tds_r))
    v_sup = sh["lines"]["meg_supply_line"]["liquid_volume_m3"]
    v_flo = sh["lines"]["flowline"]["liquid_volume_m3"]
    lean_meg, rich_meg = L.cum(lean * wm_l), L.cum(rich * wm_r)

    def held(tm):
        near = np.abs(L.t - tm) <= 60
        tk_l = (x["LT-101"][near] * sh["tanks"]["TK-L"]["floor_area_m2"]
                * R.rho(x["TT-102"][near], wm_l[near], tds_l[near]) * wm_l[near]).mean()
        tk_r = (x["LT-201"][near] * sh["tanks"]["TK-R"]["floor_area_m2"] * x["DT-401"][near] * wm_tr[near]).mean()
        back = float(np.interp(np.interp(tm, vs[1], vs[0]) - v_sup, vs[0], vs[1]))
        ahead = float(np.interp(np.interp(tm, vf[1], vf[0]) + v_flo, vf[0], vf[1]))
        return tk_l + tk_r + L.between(lean_meg, back, tm) + L.between(rich_meg, tm, ahead)

    made_up = L.between(L.cum(x["FI-601"] * purity), L.ws, L.we)
    days = (L.we - L.ws) / 1440.0
    return (made_up - (held(L.we) - held(L.ws))) / days


def readings(folder):
    out = {}
    L = R.Loop(folder)
    fw, n = R.formation_water(L)
    ml, m = R.meg_loss(L)
    out["reference"] = (fw, ml)
    out["formation water by the chloride balance only"] = (n["by_chloride"], ml)
    out["formation water by the water balance only"] = (n["by_water"], ml)
    out["formation water as the mean of the chloride and water balances"] = ((n["by_chloride"] + n["by_water"]) / 2, ml)
    out["formation water as the median of chloride, calcium and water balances"] = (
        float(np.median([n["by_chloride"], n["by_calcium"], n["by_water"]])), ml)
    L2 = R.Loop(folder)
    L2.sheet["lines"]["meg_supply_line"]["liquid_volume_m3"] = 0.0
    L2.sheet["lines"]["flowline"]["liquid_volume_m3"] = 0.0
    fw2, _ = R.formation_water(L2)
    out["no transit: the window taken at the plant"] = (fw2, ml)
    L4 = R.Loop(folder)
    L4._rich_k = 1.0
    fw4, _ = R.formation_water(L4)
    out["rich MEG meter as read"] = (fw4, ml)
    saved = R.isolated_tank_loss
    R.isolated_tank_loss = lambda L_: (0.0, {}, 0.0)
    ml3, _ = R.meg_loss(R.Loop(folder))
    R.isolated_tank_loss = saved
    out["MEG loss from the measured streams only"] = (fw, ml3)
    out["MEG loss by the inventory route (fails the examples)"] = (fw, inventory_route(R.Loop(folder)))
    return out


def main(roots):
    truth = {}
    for f in [HERE / "environment" / "data" / "expected.json", HERE / "tests" / "loop_truth.json"]:
        truth.update(json.loads(f.read_text()))
    folders = sorted(p for r in roots for p in Path(r).iterdir() if p.is_dir())
    table = {}
    for f in folders:
        t = truth[f.name]
        for name, (fw, ml) in readings(f).items():
            dev = (abs(fw / t["formation_water_t_d"] - 1) / LIM["formation_water_t_d"],
                   abs(ml / t["meg_loss_kg_d"] - 1) / LIM["meg_loss_kg_d"])
            table.setdefault(name, {})[f.name] = dev
    names = [f.name for f in folders]
    lines = ["reading | " + " | ".join(names)]
    for name, row in table.items():
        cells = []
        for n_ in names:
            a, b = row[n_]
            mark = "*" if max(a, b) > 1 else ""
            cells.append(f"{a:.1f}/{b:.1f}{mark}")
        lines.append(f"{name} | " + " | ".join(cells))
    text = "\n".join(lines)
    print(text)
    (Path(__file__).parent / "misreadings.txt").write_text(
        "worst deviation in limits, formation water / MEG loss (limits 3 % and 5 %); * = a limit exceeded\n"
        + text + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main(sys.argv[1:])
