"""Second, structurally different analysis of a campaign folder (no code shared with solution/ or sim/).

Stage 1, one scan window at a time: the window's traces are fitted with the profile-weighted mixture of
exponentials for the window's absorption of free area and free centre, its shape fixed by the current stage-2
state: every list line of the window twice, once for the gas inside the furnace (Doppler and Lorentz widths of the
reactor temperature) and once for the gas in the water-cooled arms (arm temperature), weighted by their predicted
areas. The laser profile is integrated with Gauss-Legendre nodes over +-5 FWHM, the Voigt function is the real
part of the Faddeeva function, the detector term is the single-pole response. Amplitude and offset of each trace
are free nonlinear parameters here (no projection); a free loss rate per window.
Stage 2, the whole day at once: the window areas of all runs are regressed on the zone model of the cell
(furnace length at the reactor temperature, the rest of the run's gas column in the arms at the arm temperature,
the purge gas in the arms with its own O2 fraction), free: the fill fraction (purge on), the feed mole fraction,
one reactor temperature per feed run and, without an evacuated run, the manometer zero. Stages 1 and 2 are
iterated until the stage-2 values move by less than 2e-4 (relative). Without a logged laser width the first sweep fits the width per window and their median is kept. The
wavemeter readings of each scan are replaced by their straight line in step number (the laser steps evenly).

    python second_solver.py <campaign folder> [<output json>]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.special import wofz

CLIGHT = 2.99792458e10
HC_K = 1.438776877
K_B = 1.380649e-23
MASS = 31.98983 * 1.66053907e-27
DENS_PER_TORR_K = 133.322368 / K_B / 1e6
XL, WL = np.polynomial.legendre.leggauss(48)


def read(folder):
    folder = Path(folder)
    meta = json.loads((folder / "campaign.json").read_text())
    tab = np.loadtxt(folder / "o2_aband_lines.csv", delimiter=",", skiprows=1, usecols=(3, 4, 5, 6, 7))
    qtab = np.loadtxt(folder / "o2_partition_sums.csv", delimiter=",", skiprows=1)
    rows = [ln.split(",") for ln in (folder / "ringdowns.csv").read_text().splitlines()[1:]]
    run = np.array([r[0] for r in rows])
    step = np.array([int(r[1]) for r in rows])
    nu = np.array([float(r[2]) for r in rows])
    y = np.array([[float(v) for v in r[3:]] for r in rows])
    return meta, tab, qtab, run, regrid(run, step, nu), y


def regrid(run, step, nu):
    """the laser steps evenly inside a scan: replace the wavemeter readings of each scan (consecutive step numbers,
    readings within 0.03 cm-1 of the previous one) by their least-squares line in step number"""
    out = nu.copy()
    for name in np.unique(run):
        idx = np.flatnonzero(run == name)
        idx = idx[np.argsort(step[idx])]
        start = 0
        for k in range(1, len(idx) + 1):
            if k == len(idx) or step[idx[k]] != step[idx[k - 1]] + 1 or abs(nu[idx[k]] - nu[idx[k - 1]]) > 0.03:
                seg = idx[start:k]
                if len(seg) >= 4:
                    a = np.vstack([np.ones(len(seg)), step[seg] - step[seg].mean()]).T
                    coef = np.linalg.lstsq(a, nu[seg], rcond=None)[0]
                    out[seg] = a @ coef
                start = k
    return out


def q_of(qtab, t):
    return float(np.exp(np.interp(math.log(t), np.log(qtab[:, 0]), np.log(qtab[:, 1]))))


def s_of(line, qtab, t):
    nu0, s296, e = line[0], line[1], line[2]
    return s296 * q_of(qtab, 296.0) / q_of(qtab, t) * math.exp(-HC_K * e * (1 / t - 1 / 296.0)) \
        * (1 - math.exp(-HC_K * nu0 / t)) / (1 - math.exp(-HC_K * nu0 / 296.0))


def voigt(x, sig, gam):
    z = (x + 1j * gam) / (sig * math.sqrt(2))
    return wofz(z).real / (sig * math.sqrt(2 * math.pi))


def laser(shape, fwhm):
    d = 5 * fwhm * XL
    if shape == "gaussian":
        g = np.exp(-4 * math.log(2) * (d / fwhm) ** 2)
    else:
        a = fwhm / (2 * math.acosh(math.sqrt(2)))
        g = 1 / np.cosh(d / a) ** 2
    w = WL * g
    return d, w / w.sum()


def window_fit(nu, y, t, tau, shape, blend, fw, free_width):
    """free: area of the window's blend (rows of centre, weight, sigma, gamma; weights sum to 1), centre shift,
    loss rate, and the laser FWHM if free_width"""
    n = len(nu)
    use = t >= max(3 * tau, 0.12e-6)
    tt = t[use]
    yy = y[:, use]
    ng = 4 if free_width else 3

    def model(p):
        area, k0, dc = p[0] * 1e-7, p[1] * 1e5, p[2] * 1e-3
        width = p[3] * 1e-2 if free_width else fw
        d, w = laser(shape, width)
        x = nu[:, None] + d[None, :] - dc
        a = area * sum(b[1] * voigt(x - b[0], b[2], b[3]) for b in blend)
        k = k0 + CLIGHT * a
        e = (np.exp(-k[:, :, None] * tt) - np.exp(-tt / tau)) / (1 - k[:, :, None] * tau)
        return np.tensordot(e, w, axes=([1], [0]))

    def res(p):
        m = model(p[:ng])
        amp, off = p[ng:ng + n], p[ng + n:]
        return (amp[:, None] * 1000 * m + off[:, None] * 10 - yy).ravel()

    seed = [1.0, 1.45, 0.0] + ([fw * 100] if free_width else [])
    base = model(np.array(seed))
    amp0 = np.array([np.polyfit(b, v, 1)[0] for b, v in zip(base, yy)]) / 1000
    off0 = np.array([np.polyfit(b, v, 1)[1] for b, v in zip(base, yy)]) / 10
    p0 = np.r_[seed, amp0, off0]
    sp = lil_matrix((n * len(tt), ng + 2 * n), dtype=int)
    for i in range(n):
        sp[i * len(tt):(i + 1) * len(tt), :ng] = 1
        sp[i * len(tt):(i + 1) * len(tt), ng + i] = 1
        sp[i * len(tt):(i + 1) * len(tt), ng + n + i] = 1
    r = least_squares(res, p0, jac_sparsity=sp, x_scale="jac", method="trf", max_nfev=200)
    j = r.jac.toarray() if hasattr(r.jac, "toarray") else r.jac
    s2 = 2 * r.cost / (r.fun.size - r.x.size)
    cov = np.linalg.pinv(j.T @ j) * s2
    return r.x[:ng], np.sqrt(np.diag(cov)[:ng])


def analyse(folder):
    meta, tab, qtab, run, nu, y = read(folder)
    dig = meta["digitiser"]
    t = np.arange(dig["samples_per_trace"]) * dig["sample_interval_ns"] * 1e-9
    tau = meta["detector"]["time_constant_ns"] * 1e-9
    shape = meta["laser"]["spectral_profile"]
    cav = meta["cavity"]
    length, heated = float(cav["mirror_separation_cm"]), float(cav["heated_length_cm"])
    purge = cav["mirror_purge"] != "off"
    x_purge = float(cav.get("purge_o2_mole_fraction") or 0.0) if purge else 0.0
    vac = [r["pressure_torr"] for r in meta["runs"] if r["role"] == "vacuum"]
    runs = {r["run"]: dict(r) for r in meta["runs"] if r["role"] != "vacuum" and np.any(run == r["run"])}
    samples = [k for k, v in runs.items() if v["role"] == "sample"]
    srt = np.sort(nu)
    cut = np.r_[srt[0] - 1, srt[1:][np.diff(srt) > 0.1] - 1e-6, srt[-1] + 1]
    win = np.searchsorted(cut, nu) - 1
    line_of = []
    for w in range(win.max() + 1):
        c = nu[win == w].mean()
        line_of.append(list(np.flatnonzero(np.abs(tab[:, 0] - c) < 1.0)))
    # state of stage 2
    st = dict(fill=0.85 if purge else 1.0, x=0.15, zero=(min(vac) if vac else 0.0),
              temps={k: 450.0 for k in samples})
    fwhm = meta["laser"].get("fwhm_cm")
    free_width = fwhm is None
    fwhm = 0.05 if fwhm is None else float(fwhm)

    def zone_fractions(fill):
        gas = fill * length
        hot = min(heated, gas)
        return hot / length, (gas - hot) / length, 1.0 - fill

    def gas_state(name, s):
        rr = runs[name]
        if rr["role"] == "reference":
            return rr["temperature_K"], rr["temperature_K"], rr["o2_mole_fraction"]
        return s["temps"][name], rr["arm_temperature_K"], s["x"]

    def zone_areas(name, s, w):
        """predicted integrated absorption of window w from the furnace zone and from the arm zones"""
        rr = runs[name]
        t_h, t_a, xg = gas_state(name, s)
        p = rr["pressure_torr"] - s["zero"]
        fh, fa, fp = zone_fractions(s["fill"])
        hot = fh * xg * DENS_PER_TORR_K * p / t_h * np.array([s_of(tab[li], qtab, t_h) for li in line_of[w]])
        arm = (fa * xg + fp * x_purge) * DENS_PER_TORR_K * p / t_a * \
            np.array([s_of(tab[li], qtab, t_a) for li in line_of[w]])
        return hot, arm, t_h, t_a, p

    def blend_of(name, s, w):
        hot, arm, t_h, t_a, p = zone_areas(name, s, w)
        tot = hot.sum() + arm.sum()
        out = []
        for k, li in enumerate(line_of[w]):
            for a, temp in ((hot[k], t_h), (arm[k], t_a)):
                if a > 0:
                    sig = tab[li][0] * math.sqrt(K_B * temp / MASS) / 2.99792458e8
                    gam = tab[li][3] * max(p, 0.0) / 760 * (296 / temp) ** tab[li][4]
                    out.append((tab[li][0], a / tot, sig, max(gam, 1e-12)))
        return out

    def sweep(s, width_free):
        out = {}
        for name in runs:
            for w in range(len(line_of)):
                m = (run == name) & (win == w)
                out[(name, w)] = window_fit(nu[m], y[m], t, tau, shape, blend_of(name, s, w), fwhm, width_free)
        return out

    names = (["fill"] if purge else []) + (["x"] if samples else []) + ([] if vac else ["zero"]) + \
        [f"T_{k}" for k in samples]

    def unpack(p, s):
        s = dict(s, temps=dict(s["temps"]))
        for nm, v in zip(names, p):
            if nm == "fill":
                s["fill"] = v
            elif nm == "x":
                s["x"] = v * 0.1
            elif nm == "zero":
                s["zero"] = v * 0.1
            else:
                s["temps"][nm[2:]] = v * 100
        return s

    def regress(fits, s):
        keys = list(fits)
        area = np.array([fits[k][0][0] for k in keys]) * 1e-7
        err = np.array([fits[k][1][0] for k in keys]) * 1e-7

        def res(p):
            ss = unpack(p, s)
            pred = []
            for name, w in keys:
                hot, arm, *_ = zone_areas(name, ss, w)
                pred.append(hot.sum() + arm.sum())
            return (area - np.array(pred)) / err
        p0 = [s["fill"] if nm == "fill" else s["x"] * 10 if nm == "x" else s["zero"] * 10 if nm == "zero"
              else s["temps"][nm[2:]] / 100 for nm in names]
        lo = [0.3 if nm == "fill" else 0.005 if nm == "x" else -5 if nm == "zero" else 1.5 for nm in names]
        hi = [1.0 if nm == "fill" else 10 if nm == "x" else 5 if nm == "zero" else 15 for nm in names]
        p0 = np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6)
        sol = least_squares(res, p0, bounds=(lo, hi), x_scale="jac")
        return unpack(sol.x, s)

    for it in range(8):                      # until the stage-2 state stops moving
        wf = free_width and it == 0
        fits = sweep(st, wf)
        if wf:
            fwhm = float(np.median([v[0][3] for v in fits.values()])) * 1e-2
            fits = sweep(st, False)
        new = regress(fits, st)
        moved = max([abs(new["fill"] / st["fill"] - 1), abs(new["x"] / st["x"] - 1)]
                    + [abs(new["temps"][k] / st["temps"][k] - 1) for k in samples])
        st = new
        if it >= 2 and moved < 2e-4:
            break
    return {"fill_fraction": float(st["fill"]), "o2_mole_fraction": float(st["x"]),
            "temperature_K": {k: float(st["temps"][k]) for k in samples},
            "laser_fwhm_cm": fwhm, "manometer_zero_torr": float(st["zero"])}


if __name__ == "__main__":
    out = analyse(sys.argv[1])
    txt = json.dumps(out, indent=1)
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(txt + "\n")
    print(txt)
