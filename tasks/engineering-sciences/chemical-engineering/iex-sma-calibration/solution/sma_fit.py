"""Steric-mass-action parameters of a protein on an ion-exchange column from one calibration campaign.

    python3 sma_fit.py <campaign folder> <answer.json>

Every run is rebuilt from the records: the pump program and the buffers give the salt leaving the pump; the mixer
(one well-mixed volume) and the tubing (plug flow, volume from length and bore) carry it, and the sample, to the
column; the column is the stated series of well-mixed stages with kinetic SMA binding; the outlet travels through the
stated sequence of tubes and detector cells. The UV reading is modelled as the detector forms it: protein absorbance
plus the salt-dependent baseline found in the blank run, through the stray-light response, minus the auto-zero at
the start. nu and k_eq are first fitted to the low-load gradients (which do not depend on sigma), then nu, k_eq and
sigma together to all runs, the low loads at 280 nm and the high loads at 300 nm, where the detector is still nearly
linear.
"""
import csv
import json
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

DT = 0.1          # s, grid of the extra-column calculation
CHECK_PUMP = True  # estimate the delivered share of B from the conductivity (the evidence scripts switch it off)
NOISE_mAU = 0.2


def read_campaign(folder):
    folder = Path(folder)
    system = json.loads((folder / "system.json").read_text())
    with open(folder / "runs.csv", newline="") as fh:
        runs = list(csv.DictReader(fh))
    for r in runs:
        for k in ("sample_mg_mL", "sample_NaCl_mM", "load_volume_mL", "load_flow_mL_min", "flow_mL_min", "wash_CV",
                  "start_pctB",
                  "end_pctB", "gradient_CV", "hold_CV"):
            r[k] = float(r[k])
        data = np.genfromtxt(folder / r["file"], delimiter=",", names=True)
        r["data"] = data
    std = np.genfromtxt(folder / "conductivity_standards.csv", delimiter=",", names=True)
    system["_standards"] = (std["NaCl_mM"], std["kappa_25C_mS_cm"])
    apply_log(folder, system, runs)
    return system, runs


def run_list(text, all_ids):
    """Run ids a remark names: 'R06 to R08', 'R06-R08', 'R05 onwards', or single ids."""
    ids = []
    for a, b in re.findall(r"R(\d+)\s*(?:to|-|\u2013)\s*R(\d+)", text):
        ids += [f"R{i:02d}" for i in range(int(a), int(b) + 1)]
    for a in re.findall(r"R(\d+)\s+onwards", text):
        ids += [x for x in all_ids if int(x[1:]) >= int(a)]
    if not ids:
        ids = re.findall(r"R\d+", text)
    return [x for x in ids if x in all_ids]


def apply_log(folder, system, runs):
    """Remarks of run_log.txt that change what a run was: a repacked bed, a re-measured or diluted sample."""
    path = Path(folder) / "run_log.txt"
    if not path.exists():
        return
    ids = [r["run"] for r in runs]
    by_id = {r["run"]: r for r in runs}
    A = system["buffers"]["A"]["NaCl_mM"]
    for line in path.read_text().splitlines():
        if " | " not in line:
            continue
        text = line.split(" | ", 1)[1]
        low = text.lower()
        if "repacked" in low:
            h = re.search(r"bed height ([\d.]+) cm", text)
            n = re.search(r"\((\d+) stages\)", text)
            for i in run_list(text, ids):
                col = dict(system["column"], bed_height_cm=float(h.group(1)))
                col["stages"] = int(n.group(1)) if n else int(round(system["column"]["stages"] * col["bed_height_cm"]
                                                                 / system["column"]["bed_height_cm"]))
                by_id[i]["_column"] = col
        m = re.search(r"re-measured[^:]*:\s*([\d.]+)\s*mg/mL", text)
        if m:
            for i in run_list(text, ids):
                by_id[i]["sample_mg_mL"] = float(m.group(1))
        if "diluted 1+1 with buffer a" in low:
            for i in run_list(text, ids):
                by_id[i]["sample_mg_mL"] /= 2.0
                by_id[i]["sample_NaCl_mM"] = (by_id[i]["sample_NaCl_mM"] + A) / 2.0


def measured_nacl(system, r):
    """NaCl (mM) in the conductivity cell from the monitor's reading: undo the monitor's own compensation, bring the
    conductivity to 25 C with the NaCl coefficient, then invert the standards (made up in the buffer species)."""
    d = r["data"]
    std_c, std_k = system["_standards"]
    comp = system["instrument"]["detectors"]["conductivity"]["temperature_compensation_pct_per_K"]
    T = d["cond_temperature_C"]
    kT = d["conductivity_mS_cm"] * (1 + comp / 100 * (T - 25)) if comp is not None else d["conductivity_mS_cm"]
    k25 = kT / (1 + 0.0212 * (T - 25))
    grid = np.linspace(std_c[0], std_c[-1], 20001)
    curve = PchipInterpolator(std_c, std_k)(grid)
    return np.interp(k25, curve, grid)


def tube_mL(seg):
    return np.pi * (seg["id_mm"] / 20.0) ** 2 * seg["length_cm"]


def mixed(u, V, Q, start):
    """Well-mixed volume V holding `start` at t = 0; exact for an inlet linear between grid points. Q[i] is the flow
    during the step ending at grid point i."""
    out = np.empty_like(u)
    out[0] = start
    for i in range(1, len(u)):
        a = np.exp(-Q[i] * DT / V)
        g = V / (Q[i] * DT) * (1 - a)
        u0 = u[i - 1] if i > 1 else start
        out[i] = u[i] + (out[i - 1] - u0) * a - (u[i] - u0) * g
    return out


def shift(pumped, x, V, start):
    """Plug flow through V: the value now is the one that entered V mL of pumped volume ago."""
    return np.interp(pumped - V, pumped, x, left=start)


class Run:
    """Extra-column part of one run, independent of the SMA parameters. Salt means sodium ion here: NaCl plus the
    sodium the buffer species bring."""

    def __init__(self, system, r):
        ins, buf, prot = system["instrument"], system["buffers"], system["protein"]
        mcol = system["column"]
        col = r.get("_column", mcol)
        self.r, self.prot, self.col = r, prot, col
        self.Qe = r["flow_mL_min"] / 60.0
        self.Ql = r["load_flow_mL_min"] / 60.0
        # the method converts CV with the column volume of the system file; the bed may differ after a repack
        self.Vc = np.pi * (mcol["diameter_cm"] / 2) ** 2 * mcol["bed_height_cm"]
        self.Vbed = np.pi * (col["diameter_cm"] / 2) ** 2 * col["bed_height_cm"]
        na = buf.get("Na_from_buffer_mM", 0.0)
        A, B = buf["A"]["NaCl_mM"] + na, buf["B"]["NaCl_mM"] + na
        Qe, Ql = self.Qe, self.Ql
        self.t_load = t_load = r["load_volume_mL"] / Ql
        t0g = t_load + r["wash_CV"] * self.Vc / Qe
        t_grad = r["gradient_CV"] * self.Vc / Qe
        t_end = t0g + t_grad + r["hold_CV"] * self.Vc / Qe
        t = np.arange(0.0, t_end + DT / 2, DT)
        Qarr = np.where(t <= t_load + 1e-9, Ql, Qe)
        pumped = Ql * np.minimum(t, t_load) + Qe * np.maximum(t - t_load, 0.0)
        frac = np.clip((t - t0g) / t_grad, 0.0, 1.0)
        pct = np.where(t < t0g, r["start_pctB"], r["start_pctB"] + (r["end_pctB"] - r["start_pctB"]) * frac)
        self._pump = (A, B, pct)
        c_sample = r["sample_mg_mL"] / prot["molar_mass_g_mol"] * 1000.0
        s_sample = r["sample_NaCl_mM"] + na
        load = t < t_load
        vol = {s["name"]: tube_mL(s) for s in ins["tubing"]}
        self.t, self.pumped, self.A = t, pumped, A

        def inlet(fB):
            salt = A + (B - A) * pct * fB / 100.0
            if r["sample_inlet"] == "pump_A":
                salt = np.where(load, s_sample, salt)
                ms = mixed(salt, ins["mixer_mL"], Qarr, A)
                mp = mixed(np.where(load, c_sample, 0.0), ins["mixer_mL"], Qarr, 0.0)
                vs, vp = shift(pumped, ms, vol["mixer_to_valve"], A), shift(pumped, mp, vol["mixer_to_valve"], 0.0)
            else:
                ms = mixed(salt, ins["mixer_mL"], Qarr, A)
                vs = np.where(load, s_sample, shift(pumped, ms, vol["mixer_to_valve"], A))
                vp = np.where(load, c_sample, 0.0)
            return shift(pumped, vs, vol["valve_to_column"], A), shift(pumped, vp, vol["valve_to_column"], 0.0)

        def to_cell(name):
            v = 0.0
            for item in ins["after_column"]:
                if item == name:
                    return v + ins["detectors"][name]["cell_volume_mL"] / 2
                v += vol[item] if item in vol else ins["detectors"][item]["cell_volume_mL"]
            return v

        self.V_to_uv = to_cell("UV")
        # the pump may deliver less B than set: compare the salt the conductivity monitor saw with the salt the flow
        # path would carry to its cell at the full setting (no protein), and scale the B part by least squares
        self.cin_s, self.cin_p = inlet(1.0)
        self.fB = 1.0
        if CHECK_PUMP:
            s1, _ = self.outlet(1.0, 1e-12, 0.0, system["binding"])
            at_cell = np.interp(r["data"]["time_s"], t, shift(pumped, s1, to_cell("conductivity"), A))
            seen = measured_nacl(system, r) + na
            use = at_cell - A > 0.2 * (B - A) * r["end_pctB"] / 100.0
            x, y = at_cell[use] - A, seen[use] - A
            self.fB = float(x @ y / (x @ x)) if use.any() else 1.0
            self.cin_s, self.cin_p = inlet(self.fB)
        d = r["data"]
        self.ts = d["time_s"]
        self.uv = {wl: d[f"uv{wl}_mAU"] for wl in (280, 300) if f"uv{wl}_mAU" in d.dtype.names}
        self.uvdet = ins["detectors"]["UV"]

    def outlet(self, nu, keq, sig, binding):
        col = self.col
        N = col["stages"]
        Vl = col["porosity_salt"] * self.Vbed / N
        Vlp = col["porosity_protein"] * self.Vbed / N
        Vs = (1 - col["porosity_salt"]) * self.Vbed / N
        Lam, kk = col["ionic_capacity_mM"], binding["k_kin_per_s"]
        Qe, Ql, tl = self.Qe, self.Ql, self.t_load
        t, cs_in, cp_in = self.t, self.cin_s, self.cin_p

        def f(tt, y):
            Q = Ql if tt < tl else Qe
            cs, cp, q = y[:N], y[N:2 * N], y[2 * N:]
            free = np.maximum(Lam - (nu + sig) * q, 0.0)
            rate = kk * (cp * (free / Lam) ** nu - q / keq * (np.maximum(cs, 1e-12) / Lam) ** nu)
            a = np.interp(tt, t, cs_in)
            b = np.interp(tt, t, cp_in)
            dcs = (Q * (np.concatenate(([a], cs[:-1])) - cs) + Vs * nu * rate) / Vl
            dcp = (Q * (np.concatenate(([b], cp[:-1])) - cp) - Vs * rate) / Vlp
            return np.concatenate((dcs, dcp, rate))
        y0 = np.concatenate((np.full(N, self.A), np.zeros(2 * N)))
        sol = solve_ivp(f, (0.0, t[-1]), y0, method="BDF", jac_sparsity=pattern(N), rtol=1e-6, atol=1e-10,
                        max_step=2.0, dense_output=True)
        Y = sol.sol(t)
        return Y[N - 1], Y[2 * N - 1]

    def uv_model(self, salt_out, prot_out, b, wl):
        s = shift(self.pumped, salt_out, self.V_to_uv, self.A)
        p = shift(self.pumped, prot_out, self.V_to_uv, 0.0)
        c_mg = p * self.prot["molar_mass_g_mol"] / 1000.0
        A = self.prot[f"extinction_{wl}nm_mL_mg_cm"] * c_mg * self.uvdet["path_cm"] + b[wl] * s / 1000.0
        st = self.uvdet["stray_light_fraction"]
        obs = -np.log10((10 ** (-A) + st) / (1 + st))
        return np.interp(self.ts, self.t, 1000.0 * (obs - obs[0]))


_PAT = {}


def pattern(N):
    if N not in _PAT:
        S = lil_matrix((3 * N, 3 * N))
        for i in range(N):
            for a in range(3):
                for c in range(3):
                    S[a * N + i, c * N + i] = 1
            if i:
                S[i, i - 1] = 1
                S[N + i, N + i - 1] = 1
        _PAT[N] = S.tocsr()
    return _PAT[N]


def blank_slopes(run):
    """Salt baseline of the UV from the blank gradient: reading = slope * (salt at the cell - start salt)."""
    s = shift(run.pumped, salt_without_binding(run), run.V_to_uv, run.A)
    ds = np.interp(run.ts, run.t, s - s[0]) / 1000.0  # M
    return {wl: float(np.dot(ds, run.uv[wl]) / np.dot(ds, ds)) / 1000.0 for wl in run.uv}  # AU per M


def salt_without_binding(run):
    out_s, _ = run.outlet(1.0, 1e-12, 0.0, BINDING)
    return out_s


BINDING = None
RUNS = None
SLOPES = None


def default_channel(run):
    """Low loads at 280 nm; large loads at 300 nm where the detector records it (the peaks stay well inside its
    range there), otherwise at 280 nm through the stray-light response."""
    return 300 if run.r["kind"] == "load" and 300 in run.uv else 280


CHANNEL = default_channel


def simulate(job):
    i, nu, keq, sig = job
    run = RUNS[i]
    so, po = run.outlet(nu, keq, sig, BINDING)
    wl = CHANNEL(run)
    return (run.uv_model(so, po, SLOPES, wl) - run.uv[wl]) / NOISE_mAU


def solve(folder, adjust=None, channel=None, start=None, workers=2):
    """Fit one campaign. adjust(system, runs) may change the parsed records first and channel(run) may choose the
    UV channel of a run (both used by the evidence scripts); start=(nu, k_eq, sigma) skips the search and refines
    from there (used only for the noise redraws, which measure scatter, not the search)."""
    global BINDING, RUNS, SLOPES, CHANNEL
    system, runs = read_campaign(folder)
    if adjust:
        adjust(system, runs)
    BINDING = system["binding"]
    RUNS = [Run(system, r) for r in runs]
    if channel:
        CHANNEL = channel
    blank = [r for r in RUNS if r.r["kind"] == "blank"][0]
    SLOPES = blank_slopes(blank)
    low = [i for i, r in enumerate(RUNS) if r.r["kind"] == "gradient"]
    high = [i for i, r in enumerate(RUNS) if r.r["kind"] == "load"]
    with Pool(workers) as pool:
        def resid(x, idx, sig=None):
            s = x[2] if sig is None else sig
            parts = pool.map(simulate, [(i, x[0], np.exp(x[1]), s) for i in idx])
            return np.concatenate(parts)

        if start is not None:
            f = least_squares(lambda x: resid(x, low + high), [start[0], np.log(start[1]), start[2]],
                              x_scale=[0.3, 0.3, 0.2 * start[2]], diff_step=1e-4, max_nfev=60)
            return {"nu": float(f.x[0]), "k_eq": float(np.exp(f.x[1])), "sigma": float(f.x[2]),
                    "chi2_per_point": float(2 * f.cost / len(f.fun))}
        # stage 1: nu and k_eq from the low loads, from a few starting charges
        best = None
        for nu0 in (3.5, 5.5, 7.5):
            try:
                f = least_squares(lambda x: resid(x, low, 0.0), [nu0, np.log(0.03)], x_scale=[0.5, 0.5],
                                  diff_step=1e-3, max_nfev=40)
            except Exception:  # noqa: BLE001
                continue
            if best is None or f.cost < best.cost:
                best = f
        nu1, lk1 = best.x
        # stage 2: all three, sigma started from a short scan on the largest load
        big = high[-1]
        scan = []
        for s0 in (15.0, 40.0, 100.0, 200.0):
            r = resid([nu1, lk1, s0], [big])
            scan.append((float(r @ r), s0))
        s0 = min(scan)[1]
        f = least_squares(lambda x: resid(x, low + high), [nu1, lk1, s0], x_scale=[0.3, 0.3, 0.2 * s0],
                          diff_step=1e-4, max_nfev=60)
    return {"nu": float(f.x[0]), "k_eq": float(np.exp(f.x[1])), "sigma": float(f.x[2]),
            "chi2_per_point": float(2 * f.cost / len(f.fun))}


def main():
    Path(sys.argv[2]).write_text(json.dumps(solve(sys.argv[1]), indent=1) + "\n")


if __name__ == "__main__":
    main()
