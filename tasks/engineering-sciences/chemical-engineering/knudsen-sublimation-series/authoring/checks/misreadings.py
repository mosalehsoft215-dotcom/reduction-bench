"""The reference analysis with one misreading at a time, on every series it changes.

    python3 misreadings.py > misreadings.csv          (the bundled series)
    python3 misreadings.py <folder> > out.csv         (seven series written by series.py, with truth.json)

Readings:
  none              the reference itself
  no_waits          only the programs count: the waits under vacuum before and after them left out
  waits_30_20       waits taken as 30 min before and 20 min after every program (the public series' values)
                    instead of the logged times
  spare_zero        the spare thermometer without a certificate taken as reading true
  spare_as_foil     that spare's offset taken up as a new orifice conductance from its first run on
  one_foil          one conductance for every run (a new foil without reference runs gets the old one's)
  leads_zero        a 2-wire Pt100 whose leads were not measured taken as having none
  no_transition     the sample's transition ignored (form I throughout)
  form_reverts      form II only while the crystal is above the transition temperature
  form_by_run       the form changes between runs, not where the crystal first reaches the transition
  keep_flagged      runs the log rules out kept
  shared_left_out   a stay holding two programs left out (allowed: it only loses data)
  no_dCp            the heat-capacity term left out
  cert_everywhere   a certified spare's offset applied at every temperature, also outside the range its
                    certificate covers (those stays kept)
  outside_unknown   the stays outside a certificate's range kept, the spare's error there taken as one constant
                    unknown amount (fitted)
  cert_sign         a certified spare's offset applied with the wrong sign (an error read as a correction)
  liner_dynamics    a thermometer in the liner ignored: the liner's lag and sublimation cooling applied as if the
                    block were held on the program
  check_ignored     a Pt100's check at 0 and 100 C ignored (cell.json's R0 and A used)
  check_with_leads  that check used without taking the leads off
  evap_ignored      the evaporation coefficients ignored (the crystals taken to pass vapour freely)
  evap_sample_only  only the sample's evaporation coefficient applied, not the reference material's
  weighed_part      a weighing made with a part still on the cell taken as it stands
  no_selfcooling    a crystal in a liner taken to be at the block temperature (the heat its sublimation takes
                    and the liner's heat capacity left out)
  quasi_static      the liner's sublimation cooling kept but its heat capacity (and the crystals') left out, so
                    the crystal follows the block without lag
  liner_cp_only     the liner's heat capacity kept, the crystals' own heat capacity left out
  no_radiation      the heat radiation the liner's orifice sends to the 77 K shroud left out (sublimation cooling,
                    heat capacities and lag kept)
  shroud_at_room    that radiation taken to a shroud at room temperature (295 K) instead of 77 K
  form_at_block     the sample's form decided on the block temperature instead of the crystal's
  offsets_from_hold the offset of every thermometer (spare, leads) read off the hold, where the controlled
                    thermometer shows the set temperature: none
Deviation of dH298 and p298 from the truth, in per cent and in limits (tests/effusion_truth.json).
"""
import copy
import json
import os
import sys
from pathlib import Path

import numpy as np

BUNDLE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BUNDLE / "solution"))
import analyse_series as E  # noqa: E402

limits = json.loads((BUNDLE / "tests" / "effusion_truth.json").read_text())["limits_percent"]


class Reverting(E.Model):
    """Form II only at temperatures at or above the transition."""

    def predict(self, x):
        lnC = x[2:2 + len(self.foils)]
        out = np.zeros(len(self.stays))
        for j, seg in enumerate(self.segments):
            total = 0.0
            for kind, rid, t, read, wait in seg:
                r = self.by[rid]
                T = self.kelvin(rid, read, x)
                C = np.exp(lnC[self.foils.index(r["foil"])]) * 1e-7
                if r["compound"] == "reference":
                    f = np.exp(self.ref_lnp(T)) * np.sqrt(self.Mr / (2 * np.pi * E.R * T))
                else:
                    la = E.ln_p_a(T, x[0] * 1e3, x[1], self.dCp)
                    la = np.where(T >= self.trs[0], la + E.ln_p_b_minus_a(T, *self.trs), la)
                    f = np.exp(la) * np.sqrt(self.Ms / (2 * np.pi * E.R * T))
                total += C * (f[0] * wait if kind == "wait" else np.trapezoid(f, t)) * 1e6
            out[j] = total
        return out


class ByRun(E.Model):
    """The run in which the crystal first reaches the transition counts as form I, every later run as form II."""

    def predict(self, x):
        saved = self.trs
        first = None
        for seg in self.segments:
            for kind, rid, t, read, wait in seg:
                if (kind == "program" and self.by[rid]["compound"] == "sample" and first is None
                        and self.kelvin(rid, read, x).max() >= saved[0]):
                    first = rid
        out = np.zeros(len(self.stays))
        nf = len(self.foils)
        lnC = x[2:2 + nf]
        form_b, seen = False, False
        for j, seg in enumerate(self.segments):
            total = 0.0
            for kind, rid, t, read, wait in seg:
                r = self.by[rid]
                T = self.kelvin(rid, read, x)
                C = np.exp(lnC[self.foils.index(r["foil"])]) * 1e-7
                if r["compound"] == "reference":
                    f = np.exp(self.ref_lnp(T)) * np.sqrt(self.Mr / (2 * np.pi * E.R * T))
                else:
                    la = E.ln_p_a(T, x[0] * 1e3, x[1], self.dCp)
                    if form_b:
                        la = la + E.ln_p_b_minus_a(T, *saved)
                    f = np.exp(la) * np.sqrt(self.Ms / (2 * np.pi * E.R * T))
                total += C * (f[0] * wait if kind == "wait" else np.trapezoid(f, t)) * 1e6
                if rid == first and kind == "program":
                    seen = True
            if seen:
                form_b = True
            out[j] = total
        return out


def misread(folder, reading):
    cls = {"form_reverts": Reverting, "form_by_run": ByRun}.get(reading, E.Model)
    m = cls(folder)
    if reading in ("form_reverts", "form_by_run", "no_transition") and m.trs is None:
        return None
    if reading == "no_transition":
        m.trs = None
    elif reading in ("no_waits", "waits_30_20"):
        changed = False
        for seg in m.segments:
            for k, (kind, rid, t, read, wait) in enumerate(seg):
                if kind != "wait":
                    continue
                if reading == "no_waits":
                    w = 0.0
                elif k == 0:
                    w = 1800.0
                elif k == len(seg) - 1:
                    w = 1200.0
                else:
                    w = wait
                changed |= abs(w - wait) > 1
                seg[k] = (kind, rid, t, read, w)
        if not changed:
            return None
    elif reading in ("spare_zero", "spare_as_foil"):
        if m.n_spare == 0:
            return None
        if reading == "spare_as_foil":
            for r in m.by.values():
                if r["spare"] is not None:
                    r["spare"], r["foil"] = None, r["foil"] + 100
            m.layout()
    elif reading == "one_foil":
        if len(m.foils) == 1:
            return None
        for r in m.by.values():
            r["foil"] = 0
        m.layout()
    elif reading == "leads_zero":
        if not m.lead_free:
            return None
    elif reading == "keep_flagged":
        if all(s["use"] for s in m.stays):
            return None
        for s in m.stays:
            s["use"] = True
        m.layout()
    elif reading == "shared_left_out":
        if all(len(s["runs"]) == 1 for s in m.stays):
            return None
        for s in m.stays:
            if len(s["runs"]) > 1:
                s["use"] = False
        m.layout()
    elif reading == "no_dCp":
        m.dCp = 0.0
    elif reading == "cert_sign":
        if not any(r["cert"] for r in m.by.values()):
            return None
        for r in m.by.values():
            r["known_offset"] = -r["known_offset"]
        m.layout()
    elif reading in ("cert_everywhere", "outside_unknown"):
        flagged = [s_ for s_ in m.stays if not s_["use"] and all(m.by[i]["use"] for i in s_["runs"])]
        if not flagged:
            return None
        for s_ in flagged:
            s_["use"] = True
            for i in s_["runs"]:
                m.by[i]["cert"] = None
                if reading == "outside_unknown":
                    m.by[i]["spare"] = m.n_spare
        m.layout()
    elif reading == "liner_dynamics":
        if m.th.get("location") != "liner" or "liner" not in m.cell:
            return None
        m.G, m.Cl = m.cell["liner"]["conductance_W_per_K"], m.cell["liner"]["heat_capacity_J_per_K"]
    elif reading in ("evap_ignored", "evap_sample_only"):
        if all(v is None for v in m.surf.values()):
            return None
        m.surf["reference"] = None
        if reading == "evap_ignored":
            m.surf["sample"] = None
    elif reading in ("check_ignored", "check_with_leads"):
        if not m.th.get("check_ohm"):
            return None
        if reading == "check_ignored":
            m.th["check_ohm"] = None
        else:
            r1, r2 = m.th["check_ohm"]
            m.th["check_ohm"] = (r1 + m.lead, r2 + m.lead)
    elif reading == "weighed_part":
        if not any(s_.get("extra_mg") for s_ in m.stays):
            return None
        for s_ in m.stays:
            s_["after"] = s_["after"] + s_.get("extra_mg", 0.0)
        m.layout()
    elif reading == "no_selfcooling":
        if m.G is None:
            return None
        m.G, m.Cl = None, None
    elif reading == "quasi_static":
        if m.Cl is None:
            return None
        m.Cl = None
    elif reading == "liner_cp_only":
        if m.Cl is None:
            return None
        m.Cp_ref, m.Cp_smp = 0.0, 0.0
    elif reading == "offsets_from_hold":
        if m.n_spare == 0 and not m.lead_free:
            return None
    elif reading in ("no_radiation", "shroud_at_room"):
        if m.Cl is None:
            return None
        if reading == "no_radiation":
            m.rad = 0.0
        else:
            m.T_shroud = 295.0
    elif reading == "form_at_block":
        if m.Cl is None or m.trs is None:
            return None
        m.form_from_block = True
    fixed = {2 + len(m.foils) + k: 0.0 for k in range(m.n_spare)} if reading == "spare_zero" else None
    if reading == "leads_zero":
        fixed = {2 + len(m.foils) + m.n_spare: 0.0}
    if reading == "offsets_from_hold":
        fixed = {2 + len(m.foils) + k: 0.0 for k in range(m.n_spare + m.lead_free)}
    x, chi2 = E.fit(m, fixed)
    return float(x[0]), float(np.exp(x[1])), chi2


READINGS = ["none", "no_waits", "waits_30_20", "spare_zero", "spare_as_foil", "one_foil", "leads_zero",
            "no_transition", "form_reverts", "form_by_run", "keep_flagged", "shared_left_out", "no_dCp",
            "no_selfcooling", "offsets_from_hold", "cert_everywhere", "outside_unknown", "cert_sign",
            "liner_dynamics", "check_ignored", "check_with_leads", "evap_ignored", "evap_sample_only",
            "weighed_part", "quasi_static",
            "liner_cp_only", "no_radiation", "shroud_at_room", "form_at_block"]


def main():
    if len(sys.argv) > 1:
        base = Path(sys.argv[1])
        folders = {p.name: p for p in base.iterdir() if p.is_dir()}
        truth = json.loads((base / "truth.json").read_text())
    else:
        folders = {p.name: p for d in (BUNDLE / "environment/data/series", BUNDLE / "tests/series")
                   for p in d.iterdir() if p.is_dir()}
        truth = json.loads((BUNDLE / "authoring/sim/truth.json").read_text())
    only = [x for x in os.environ.get("MIS_ONLY", "").split(",") if x]
    if not os.environ.get("MIS_NOHEAD"):
        print("reading,series,dH_dev_pct,p_dev_pct,dH_in_limits,p_in_limits,worst_in_limits,chi2_per_stay")
    for reading in (only or READINGS):
        for name, t in truth.items():
            got = misread(folders[name], reading)
            if got is None:
                continue
            dH, p, chi2 = got
            a = 100 * (dH / t["dH_sub_298_kJ_mol"] - 1)
            b = 100 * (p / t["p_298_Pa"] - 1)
            la, lb = abs(a) / limits["dH_sub_298_kJ_mol"], abs(b) / limits["p_298_Pa"]
            print(f"{reading},{name},{a:+.3f},{b:+.3f},{la:.2f},{lb:.2f},{max(la, lb):.2f},{chi2:.2f}", flush=True)


if __name__ == "__main__":
    main()
