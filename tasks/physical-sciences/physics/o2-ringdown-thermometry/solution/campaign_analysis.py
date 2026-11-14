"""O2 A-band pulsed cavity ring-down: reactor temperatures, O2 mole fraction and fill fraction of one campaign
folder.

    python campaign_analysis.py <campaign folder> <output json>

Model of every averaged trace: the laser pulse excites all frequencies of its spectral profile at once; each
frequency nu' decays with its own rate k0 + c * a(nu'), so the cavity output is the profile-weighted sum of
exponentials, seen through the single-pole detector:

    y(t) = A * sum_j w_j (exp(-k_j t) - exp(-t / tau_d)) / (1 - k_j tau_d) + B

The absorption a(nu) is a sum over the zones of the cell, each a fraction of the mirror separation L:
  * the run's gas inside the furnace (heated length Lh) at the reactor temperature,
  * the run's gas that reaches past the furnace into the water-cooled arms, fill * L - Lh, at the arm temperature,
  * the purge gas, (1 - fill) * L, at the arm temperature (it absorbs only when the purge gas holds O2).
A reference run has the furnace off, so every zone is at the reference temperature. The gas in the cooled arms is
a few percent of the column but much denser and colder than the reactor gas: left out, it pulls the Boltzmann
temperature of a hot run down by several percent without spoiling the fit.

Steps:
 1. pressures: the manometer reading of the evacuated run is its zero offset and is subtracted from every run; a
    day without an evacuated run has the offset as a free parameter, fixed by runs of known gas at different fill
    pressures;
 2. windows: steps are grouped into scan windows; each window gets the lines of the list within 1 cm-1; the laser
    frequency of each step comes from a straight line through the wavemeter readings of its scan (the
    readings scatter by about 0.001 cm-1 about the even laser steps);
 3. start values from a quick pass (single exponentials, rate areas, Boltzmann slope);
 4. one joint fit of all traces of the gas runs: fill fraction shared (1 with the purge off), the O2 mole fraction
    shared by the sample runs, one reactor temperature per sample run, one loss rate per run and window, the laser
    width when the campaign has none logged, the manometer zero when there is no evacuated run; amplitude and offset
    of every trace are solved linearly inside the residual (variable projection).
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.special import voigt_profile

C_CM = 2.99792458e10
C2 = 1.438776877
KB = 1.380649e-23
M_O2 = 31.98983 * 1.66053907e-27
N_PER_TORR = 133.322368 / KB * 1e-6


class Campaign:
    def __init__(self, folder):
        folder = Path(folder)
        self.meta = json.loads((folder / "campaign.json").read_text())
        lines = np.genfromtxt(folder / "o2_aband_lines.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
        self.lines = np.column_stack([lines["wavenumber_cm"], lines["intensity_296K_cm_per_molecule"],
                                      lines["lower_energy_cm"], lines["gamma_air_cm_per_atm"], lines["n_air"]])
        qt = np.loadtxt(folder / "o2_partition_sums.csv", delimiter=",", skiprows=1)
        self.q = CubicSpline(qt[:, 0], qt[:, 1])
        raw = np.genfromtxt(folder / "ringdowns.csv", delimiter=",", skip_header=1, dtype=str)
        self.run_of = raw[:, 0]
        self.step = raw[:, 1].astype(float)
        self.nu_read = raw[:, 2].astype(float)
        self.y = raw[:, 3:].astype(float)
        dig = self.meta["digitiser"]
        self.t = np.arange(dig["samples_per_trace"]) * dig["sample_interval_ns"] * 1e-9
        self.tau = self.meta["detector"]["time_constant_ns"] * 1e-9
        self.keep = self.t >= max(3 * self.tau, 0.12e-6)
        self.tk = self.t[self.keep]
        self.runs = [r for r in self.meta["runs"] if r["role"] == "vacuum" or np.any(self.run_of == r["run"])]
        zero = [r["pressure_torr"] for r in self.runs if r["role"] == "vacuum"]
        self.gauge_zero = float(np.median(zero)) if zero else None     # None: fitted
        cav = self.meta["cavity"]
        self.length = float(cav["mirror_separation_cm"])
        self.heated = float(cav["heated_length_cm"])
        self.purge = cav["mirror_purge"] != "off"
        self.purge_x = float(cav.get("purge_o2_mole_fraction") or 0.0) if self.purge else 0.0
        fw = self.meta["laser"].get("fwhm_cm")
        self.fwhm = float(fw) if fw is not None else None              # None: fitted
        self.shape = self.meta["laser"]["spectral_profile"]
        # windows
        order = np.argsort(self.nu_read)
        gaps = np.diff(self.nu_read[order]) > 0.1
        lab = np.concatenate([[0], np.cumsum(gaps)])
        self.window = np.empty(len(self.nu_read), int)
        self.window[order] = lab
        self.nwin = lab.max() + 1
        # the laser is stepped on an even grid within each window; the wavemeter readings scatter about it, so the
        # frequency of every step is taken from a straight line through the readings of its scan
        self.nu = self.nu_read.copy()
        for run in np.unique(self.run_of):
            idx = np.flatnonzero(self.run_of == run)
            idx = idx[np.argsort(self.step[idx])]
            cut = np.flatnonzero((np.abs(np.diff(self.nu_read[idx])) > 0.05) | (np.diff(self.step[idx]) != 1)) + 1
            for seg in np.split(idx, cut):          # one contiguous scan
                if len(seg) >= 4:
                    c = np.polyfit(self.step[seg], self.nu_read[seg], 1)
                    self.nu[seg] = np.polyval(c, self.step[seg])
        self.win_lines = []
        for w in range(self.nwin):
            c = self.nu[self.window == w].mean()
            self.win_lines.append(np.flatnonzero(np.abs(self.lines[:, 0] - c) < 1.0))

    def intensity(self, t):
        l = self.lines
        return l[:, 1] * self.q(296.0) / self.q(t) * np.exp(-C2 * l[:, 2] * (1 / t - 1 / 296.0)) \
            * (1 - np.exp(-C2 * l[:, 0] / t)) / (1 - np.exp(-C2 * l[:, 0] / 296.0))

    def alpha(self, nu, w, t, p, x):
        idx = self.win_lines[w]
        l = self.lines[idx]
        s = self.intensity(t)[idx]
        sd = l[:, 0] * np.sqrt(KB * t / M_O2) / 2.99792458e8
        gl = l[:, 3] * (max(p, 0.0) / 760.0) * (296.0 / t) ** l[:, 4]
        out = np.zeros_like(nu)
        for i in range(len(idx)):
            out += s[i] * voigt_profile(nu - l[i, 0], sd[i], max(gl[i], 1e-12))
        return N_PER_TORR * p / t * x * out

    def zones(self, fill):
        """fractions of L: run gas inside the furnace, run gas in the cooled arms, purge gas"""
        gas = fill * self.length
        hot = min(self.heated, gas)
        return hot / self.length, (gas - hot) / self.length, 1.0 - fill


def laser_nodes(shape, fwhm):
    if shape == "gaussian":
        u, w = np.polynomial.hermite_e.hermegauss(41)
        return u * fwhm / (2 * np.sqrt(2 * np.log(2))), w / w.sum()
    if shape == "sech2":
        a = fwhm / (2 * np.arccosh(np.sqrt(2)))
        u = np.linspace(-9 * a, 9 * a, 181)
        w = 1 / np.cosh(u / a) ** 2
        return u, w / w.sum()
    raise ValueError(f"unknown laser profile {shape}")


def project(y, s):
    sm, ym = s.mean(1, keepdims=True), y.mean(1, keepdims=True)
    a = np.sum((s - sm) * (y - ym), 1, keepdims=True) / np.sum((s - sm) ** 2, 1, keepdims=True)
    return y - a * s - (ym - a * sm)


def single_rate(t, y):
    def f(q):
        return q[0] * np.exp(-q[1] * 1e5 * t) + q[2] - y
    q0 = [y[0] - y[-20:].mean(), 1.5, y[-20:].mean()]
    return least_squares(f, q0, x_scale=[100, 0.1, 1]).x[1] * 1e5


class JointFit:
    # global parameters and their scales; the ones a campaign does not need are left out
    SCALE = {"fwhm": 0.01, "fill": 1.0, "x": 0.1, "p0": 0.1}

    def __init__(self, cam):
        self.c = cam
        self.samples = [r["run"] for r in cam.runs if r["role"] == "sample"]
        self.ref = [r for r in cam.runs if r["role"] == "reference"]
        self.groups = []                 # (run dict, window, row indices)
        for r in cam.runs:
            if r["role"] == "vacuum":
                continue
            for w in range(cam.nwin):
                rows = np.flatnonzero((cam.run_of == r["run"]) & (cam.window == w))
                if len(rows):
                    self.groups.append((r, w, rows))
        self.y = cam.y[:, cam.keep]
        self.ng = len(self.groups)
        names = []
        if cam.fwhm is None:
            names.append("fwhm")
        if cam.purge:
            names.append("fill")
        if self.samples:
            names.append("x")
        if cam.gauge_zero is None:
            names.append("p0")
        names += [f"T_{s}" for s in self.samples]
        self.names = names
        self.nglob = len(names)
        self._arm_cache = {}

    def unpack(self, q):
        v = dict(zip(self.names, q[:self.nglob]))
        cam = self.c
        fwhm = v["fwhm"] * 0.01 if "fwhm" in v else cam.fwhm
        fill = v["fill"] if "fill" in v else 1.0
        x = v.get("x", 1.0) * 0.1
        p0 = v["p0"] * 0.1 if "p0" in v else cam.gauge_zero
        temps = {s: v[f"T_{s}"] * 100 for s in self.samples}
        k0 = q[self.nglob:] * 1e5
        return dict(fwhm=fwhm, fill=fill, x=x, p0=p0, temps=temps, k0=k0)

    def gas(self, r, par):
        """reactor-zone temperature, arm temperature and O2 fraction of the run's gas"""
        t_arm = r.get("arm_temperature_K") or r.get("temperature_K") or 296.0
        if r["role"] == "reference":
            return r["temperature_K"], r["temperature_K"], r["o2_mole_fraction"]
        return par["temps"][r["run"]], t_arm, par["x"]

    def absorption(self, g, nus, par):
        cam = self.c
        r, w, rows = self.groups[g]
        t_hot, t_arm, x_gas = self.gas(r, par)
        p = r["pressure_torr"] - par["p0"]
        f_hot, f_arm, f_purge = cam.zones(par["fill"])
        a_hot = cam.alpha(nus.ravel(), w, t_hot, p, 1.0).reshape(nus.shape)
        key = (g, round(p, 12), round(par["fwhm"], 12))
        a_arm = self._arm_cache.get(key)
        if a_arm is None:
            a_arm = cam.alpha(nus.ravel(), w, t_arm, p, 1.0).reshape(nus.shape)
            if len(self._arm_cache) > 4 * self.ng:
                self._arm_cache.clear()
            self._arm_cache[key] = a_arm
        return f_hot * x_gas * a_hot + (f_arm * x_gas + f_purge * cam.purge_x) * a_arm

    def residual(self, q):
        cam = self.c
        par = self.unpack(q)
        u, wt = laser_nodes(cam.shape, par["fwhm"])
        out = []
        for g, (r, w, rows) in enumerate(self.groups):
            nus = cam.nu[rows][:, None] + u[None, :]
            k = par["k0"][g] + C_CM * self.absorption(g, nus, par)
            e = (np.exp(-k[:, :, None] * cam.tk[None, None, :]) - np.exp(-cam.tk / cam.tau)[None, None, :]) \
                / (1 - k[:, :, None] * cam.tau)
            s = np.einsum("j,ijt->it", wt, e)
            out.append(project(self.y[rows], s).ravel())
        return np.concatenate(out)

    def sparsity(self):
        sizes = [len(rows) * self.y.shape[1] for _, _, rows in self.groups]
        m = lil_matrix((sum(sizes), self.nglob + self.ng), dtype=int)
        start = 0
        for g, size in enumerate(sizes):
            r = self.groups[g][0]
            for j, name in enumerate(self.names):
                if name in ("fwhm", "fill", "p0") or (name == "x" and r["role"] == "sample") \
                        or name == f"T_{r['run']}":
                    m[start:start + size, j] = 1
            m[start:start + size, self.nglob + g] = 1
            start += size
        return m

    def starts(self):
        """quick pass: single-exponential rates, loss from the window edges, rate areas -> Boltzmann"""
        cam = self.c
        k0, areas = [], {}
        for g, (r, w, rows) in enumerate(self.groups):
            kk = np.array([single_rate(cam.tk, yy) for yy in self.y[rows]])
            nu = cam.nu[rows]
            o = np.argsort(nu)
            base = np.median(np.r_[kk[o][:3], kk[o][-3:]])
            k0.append(base)
            areas[(r["run"], w)] = np.trapezoid(kk[o] - base, nu[o]) / C_CM
        k0 = np.array(k0)
        p0 = cam.gauge_zero
        if p0 is None:
            p0 = 0.0
            same = {}
            for r in self.ref:
                same.setdefault((r["o2_mole_fraction"], r["temperature_K"]), []).append(r)
            for rs in same.values():
                if len(rs) >= 2:
                    rs = sorted(rs, key=lambda r: r["pressure_torr"])
                    a = [np.sum([areas[(r["run"], w)] for w in range(cam.nwin)]) for r in (rs[0], rs[-1])]
                    p_lo, p_hi = rs[0]["pressure_torr"], rs[-1]["pressure_torr"]
                    if a[1] != a[0]:
                        p0 = float(np.clip((a[1] * p_lo - a[0] * p_hi) / (a[1] - a[0]), -0.5, 0.5))
                    break
        pres = lambda r: r["pressure_torr"] - p0
        fill = 1.0
        if cam.purge:
            fill = 0.9
            if self.ref:
                r = self.ref[0]
                pred = np.array([cam.intensity(r["temperature_K"])[cam.win_lines[w]].sum() for w in range(cam.nwin)])
                meas = np.array([areas[(r["run"], w)] for w in range(cam.nwin)])
                ratio = np.median(meas / pred) / (N_PER_TORR * pres(r) / r["temperature_K"])
                if abs(r["o2_mole_fraction"] - cam.purge_x) > 1e-3:
                    fill = (ratio - cam.purge_x) / (r["o2_mole_fraction"] - cam.purge_x)
            fill = float(np.clip(fill, 0.5, 1.0))
        temps, xs = [], []
        for name in self.samples:
            r = next(rr for rr in cam.runs if rr["run"] == name)
            best = None
            for t_try in np.arange(250.0, 1000.0, 5.0):
                pred = np.array([cam.intensity(t_try)[cam.win_lines[w]].sum() for w in range(cam.nwin)])
                meas = np.array([areas[(name, w)] for w in range(cam.nwin)])
                ratio = meas / pred
                spread = np.std(np.log(np.abs(ratio) + 1e-30))
                if best is None or spread < best[0]:
                    best = (spread, t_try, np.median(ratio))
            temps.append(best[1])
            xs.append(best[2] / (N_PER_TORR * pres(r) / best[1]))
        x = float(np.clip(np.median(xs) / fill, 0.005, 0.9)) if xs else 0.1
        glob = {"fwhm": 5.0, "fill": fill, "x": 10 * x, "p0": 10 * p0}
        glob.update({f"T_{s}": t / 100 for s, t in zip(self.samples, temps)})
        return np.r_[[glob[n] for n in self.names], k0 / 1e5]

    def bounds(self):
        lo = {"fwhm": 1.5, "fill": 0.3, "x": 0.005, "p0": -5.0}
        hi = {"fwhm": 15.0, "fill": 1.0, "x": 10.0, "p0": 5.0}
        lo_v = [lo.get(n, 1.5) for n in self.names] + [0.01] * self.ng
        hi_v = [hi.get(n, 15.0) for n in self.names] + [50.0] * self.ng
        return np.array(lo_v), np.array(hi_v)

    def fit(self):
        q0 = self.starts()
        lo, hi = self.bounds()
        sp = self.sparsity()
        tries = [q0]
        if "fwhm" in self.names:                     # laser width not logged: a few starting widths
            j = self.names.index("fwhm")
            tries = []
            for fw in (3.5, 5.0, 7.0):
                q = q0.copy()
                q[j] = fw
                tries.append(q)
        best = None
        for q in tries:
            q = np.clip(q, lo + 1e-6, hi - 1e-6)
            f = least_squares(self.residual, q, bounds=(lo, hi), jac_sparsity=sp, x_scale="jac",
                              diff_step=1e-4, max_nfev=60)
            if best is None or f.cost < best.cost:
                best = f
        f = least_squares(self.residual, best.x, bounds=(lo, hi), jac_sparsity=sp, x_scale="jac",
                          diff_step=1e-5, max_nfev=200)
        sigma2 = 2 * f.cost / (f.fun.size - f.x.size)
        return f.x, sigma2


def analyse(folder):
    cam = Campaign(folder)
    jf = JointFit(cam)
    q, s2 = jf.fit()
    par = jf.unpack(q)
    return {"fill_fraction": round(float(par["fill"]), 6),
            "o2_mole_fraction": round(float(par["x"]), 6),
            "temperature_K": {k: round(float(v), 3) for k, v in par["temps"].items()},
            "laser_fwhm_cm": round(float(par["fwhm"]), 6),
            "manometer_zero_torr": round(float(par["p0"]), 5),
            "noise_counts": round(float(np.sqrt(s2)), 3)}


if __name__ == "__main__":
    res = analyse(sys.argv[1])
    Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
    Path(sys.argv[2]).write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res))
