"""Reduces one Devanathan-Stachurski permeation test to D at 25 C and the sub-surface concentration C0.

    python3 permeation_reduce.py <membrane directory> <result.json>

Model: Fickian diffusion through thickness L (mean of the micrometer readings), entry face at C0 (charging on) or 0
(off), exit face with a finite oxidation rate, -D dC/dx = k C (q = D / (k L), fitted separately at each temperature; 0 for a perfect exit), solved exactly by
eigenfunction series (beta cos beta * q + sin beta = 0).  The first charging (trap filling) and the temperature change
are skipped; decay 1 and rise 2 are at the first temperature, decay 2 and rise 3 at the second (thermocouple means).
Gas bubbles that cling to the charged face for a while (the entry potential steps further past its steady charging
value while they sit there) cut the entry concentration over those stretches; the response of the exit flux to each
such stretch is added as a pair of opposite lagged steps with its own fitted loss, so its after-effect at the exit is
modelled rather than left in the transients.
One weighted least-squares fit over those four transients and the opening background gives D(T1), D(T2), C0, q(T1), q(T2)
and a linear background.  The activation energy follows from D(T1), D(T2) and the two temperatures, and gives D at 25 C.
"""
import csv
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import brentq, least_squares

F, R_GAS = 96485.33212, 8.314462618
N_TERMS = 80


def roots(q):
    """First N_TERMS positive roots of q beta cos(beta) + sin(beta) = 0 (kept for the last few q asked for)."""
    key = (float(q), N_TERMS)
    if key not in _ROOTS:
        if len(_ROOTS) > 64:
            _ROOTS.clear()
        _ROOTS[key] = _roots(float(q))
    return _ROOTS[key]


_ROOTS = {}


def _roots(q):
    if q <= 1e-9:
        return np.pi * np.arange(1, N_TERMS + 1)
    out = []
    for n in range(1, N_TERMS + 1):
        lo, hi = (n - 0.5) * np.pi, n * np.pi
        out.append(brentq(lambda b: q * b * math.cos(b) + math.sin(b), lo, hi))
    return np.array(out)


def linear_coeffs(beta, a0, a1):
    """Coefficients of A + B x/L (x/L in [0, 1]) on sin(beta x / L), with norms."""
    s, c = np.sin(beta), np.cos(beta)
    norm = 0.5 * (1 - np.sin(2 * beta) / (2 * beta))
    int_one = (1 - c) / beta
    int_x = (s - beta * c) / beta ** 2
    return (a0 * int_one + a1 * int_x) / norm


def pair_flux(t_dec, t_rise, t_on, d, phi, c0, q, L, tau=0.0):
    """Exit flux (mol m^-2 s^-1) for a decay starting at 0 from the 'on' steady state, then a rise starting at t_on.
    d is the lattice diffusivity, phi the reversible-trap capacity (transients run at d / phi), tau the entry-surface
    time constant: the boundary follows the charging through a first-order lag, so every mode e^(-lam t) of the step
    response becomes (e^(-lam t) - e^(-t/tau)) / (1 - lam tau) and the steady change is reached as 1 - e^(-t/tau)."""
    beta = roots(q)
    lam = d / phi * beta ** 2 / L ** 2
    slope_on = -c0 / (1 + q)                    # steady 'on' profile: C0 + slope_on * x/L
    a = linear_coeffs(beta, c0, slope_on)       # decay: u = C (steady 'off' is zero)
    j_ss = -(d / L) * slope_on
    a_end = a * np.exp(-lam * t_on)             # profile at the switch, relative to zero
    b = a_end - linear_coeffs(beta, c0, slope_on)   # rise: u = C - steady 'on'
    g_dec = -(d / L) * a * beta * np.cos(beta)
    g_rise = -(d / L) * b * beta * np.cos(beta)
    tau = max(tau, 1e-3)
    den = 1 - lam * tau
    den = np.where(np.abs(den) < 1e-9, 1e-9, den)

    def lagged(tt, j0, jinf, g):
        lag = np.exp(-tt / tau)
        modes = (np.exp(-np.outer(tt, lam)) - lag[:, None]) / den[None, :]
        return j0 + (jinf - j0) * (1 - lag) + modes @ g

    j_dec = lagged(t_dec, j_ss, 0.0, g_dec)
    j_on = float(g_dec @ np.exp(-lam * t_on))   # flux at the switch (the series at t = 0 itself does not converge)
    j_rise = lagged(t_rise, j_on, j_ss, g_rise)
    return j_dec, j_rise


def step_flux(tt, d, phi, c0, q, L, tau=0.0):
    """Exit flux for the entry target stepped from zero to c0 at tt = 0 (through the entry lag), from an empty membrane;
    zero before the step."""
    beta = roots(q)
    lam = d / phi * beta ** 2 / L ** 2
    slope_on = -c0 / (1 + q)
    g = (d / L) * linear_coeffs(beta, c0, slope_on) * beta * np.cos(beta)
    j_ss = -(d / L) * slope_on
    tau = max(tau, 1e-3)
    den = 1 - lam * tau
    den = np.where(np.abs(den) < 1e-9, 1e-9, den)
    x = np.clip(tt, 0.0, None)
    lag = np.exp(-x / tau)
    modes = (np.exp(-np.outer(x, lam)) - lag[:, None]) / den[None, :]
    return np.where(tt > 0, j_ss * (1 - lag) + modes @ g, 0.0)


def bubble_windows(t, pot, stretches, tau):
    """Stretches of a charging period where the entry potential sits beyond its steady charging value (by more than
    0.9 mV in the direction charging moves it, on a 15-reading running median, for 30 s or more): bubbles shading part
    of the face.  Returns (start, end) pairs."""
    smooth = median_filter(pot, size=15, mode="nearest")
    off_level = np.median(pot[t < stretches[0][0]]) if (t < stretches[0][0]).any() else pot[0]
    out = []
    for a, b in stretches:
        sel = np.flatnonzero((t > a + 5 * tau + 60.0) & (t < b))
        if sel.size < 30:
            continue
        on_level = np.median(pot[sel])
        sign = 1.0 if on_level >= off_level else -1.0
        beyond = (smooth[sel] - on_level) * sign > 0.9
        edges = np.flatnonzero(np.diff(np.r_[0, beyond.astype(int), 0]))
        for i0, i1 in zip(edges[::2], edges[1::2]):
            if i1 - i0 >= 15:
                dt = t[sel[1]] - t[sel[0]]
                out.append((t[sel[i0]] - dt / 2, t[sel[i1 - 1]] + dt / 2))
    return out


def entry_time_constant(t, pot, switches, skip=()):
    """The entry potential follows the sub-surface concentration: after each switch it relaxes exponentially with the
    entry surface's time constant. One exponential fit per switch (over up to 3000 s), median of the four."""
    taus = []
    for k, sw in enumerate(switches):
        nxt = switches[k + 1] if k + 1 < len(switches) else t[-1]
        sel = (t > sw) & (t < min(nxt, sw + 3000.0))
        for a, b in skip:
            sel &= ~((t >= a) & (t < b))
        x, y = t[sel] - sw, pot[sel]
        if x.size < 20:
            continue
        before = pot[(t < sw) & (t > sw - 60.0)].mean()
        if abs(y[0] - before) > 0.9 * abs(y[-20:].mean() - before):
            taus.append(1e-3)                          # already settled at the first reading: a clean surface
            continue
        def res(p):
            return p[0] + (p[1] - p[0]) * np.exp(-x / np.exp(p[2])) - y
        sol = least_squares(res, [y[-20:].mean(), y[:2].mean(), math.log(20.0)],
                            bounds=([-np.inf, -np.inf, math.log(1e-3)], [np.inf, np.inf, math.log(5000.0)]))
        taus.append(math.exp(sol.x[2]))
    return float(np.median(taus)) if taus else 1e-3


def read(folder, as_cut=False):
    with open(folder / "transient.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["time_s"]) for r in rows])
    cur = np.array([float(r["exit_current_ua"]) for r in rows])
    temp = np.array([float(r["cell_thermocouple_c"]) for r in rows])
    pot = np.array([float(r["entry_potential_mv"]) for r in rows])
    with open(folder / "events.csv", newline="") as fh:
        ev = [(float(r["time_s"]), r["event"]) for r in csv.DictReader(fh)]
    offs = [e[0] for e in ev if e[1].startswith("charging off")]
    ons = [e[0] for e in ev if e[1].startswith("charging on")]
    with open(folder / "micrometer.csv", newline="") as fh:
        L = np.mean([float(r["thickness_mm"]) for r in csv.DictReader(fh)]) * 1e-3
    sheet = json.loads((folder / "specimen.json").read_text())
    if as_cut:
        L = sheet["thickness_as_cut_mm"] * 1e-3
    return t, cur, temp, offs, ons, L, sheet["exposed_area_cm2"] * 1e-4, sheet, pot


N_LATTICE = 8.5e5


def reduce_membrane(folder, exit_fit=True, use_setpoints=False, as_cut=False, ea_fixed=None, entry=True, screen=True,
                    c0_split=False, traps="fit", bubbles=True):
    """The keyword options switch single steps to common shortcuts (authoring/evidence/bench.py uses them)."""
    t, cur, temp, offs, ons, L, area, sheet, pot = read(Path(folder), as_cut)
    switches = [offs[0], ons[1], offs[1], ons[2]]
    tau_entry = entry_time_constant(t, pot, switches) if entry else 0.0
    with open(Path(folder) / "events.csv", newline="") as fh:
        hold = [float(r["time_s"]) for r in csv.DictReader(fh) if r["event"].startswith("thermostat")][0]
    shades = []
    if bubbles:
        shades = bubble_windows(t, pot, [(ons[1], hold), (hold, offs[1]), (ons[2], t[-1] + 1)], max(tau_entry, 1.0))
        if shades and entry:
            tau_entry = entry_time_constant(t, pot, switches, shades)
    stride = max(1, t.size // 3500)                     # about 3500 readings are plenty for the fit
    t, cur, temp = t[::stride], cur[::stride], temp[::stride]
    n_sheet = sheet.get("reversible_trap_density_mol_per_m3", 0.0) if traps else 0.0
    e_b = sheet.get("trap_binding_energy_kj_per_mol", 0.0) * 1e3
    free_traps = traps == "fit" and e_b > 0                 # trap density from the two temperatures (E_b of the sheet)
    capacity = lambda T, n=n_sheet: 1 + n / N_LATTICE * math.exp(e_b / (R_GAS * T))
    end = t[-1] + 1
    win = {"bg": t < ons[0], "d1": (t >= offs[0]) & (t < ons[1]), "r2": (t >= ons[1]) & (t < offs[1] - 1e-6),
           "d2": (t >= offs[1]) & (t < ons[2]), "r3": t >= ons[2]}
    # rise 2 ends where the thermostat change starts: the event list has it between rise 2 and decay 2
    win["r2"] = (t >= ons[1]) & (t < hold)
    T1 = temp[win["d1"] | win["r2"]].mean() + 273.15
    T2 = temp[win["d2"] | win["r3"]].mean() + 273.15
    if use_setpoints:
        T1, T2 = (v + 273.15 for v in json.loads((Path(folder) / "specimen.json").read_text())["thermostat_setpoints_c"])
    scale = F * area * 1e6                      # mol m^-2 s^-1 -> uA

    cap1, cap2 = capacity(T1), capacity(T2)
    phi1, phi2 = cap1, cap2

    def model(p):
        ld1, ld2, lc0, lc0b, q1, q2, b0, b1, lnt = p[:9]
        d1, d2, c0 = np.exp(ld1), np.exp(ld2), np.exp(lc0)
        phi1, phi2 = (capacity(T1, math.exp(lnt)), capacity(T2, math.exp(lnt))) if free_traps else (cap1, cap2)
        c0b = np.exp(lc0b) if c0_split else c0
        tau = tau_entry
        pred = b0 + b1 * t / 3600
        out = np.full(t.size, np.nan)
        j1d, j1r = pair_flux(t[win["d1"]] - offs[0], t[win["r2"]] - ons[1], ons[1] - offs[0], d1, phi1, c0, q1, L, tau)
        j2d, j2r = pair_flux(t[win["d2"]] - offs[1], t[win["r3"]] - ons[2], ons[2] - offs[1], d2, phi2, c0b, q2, L,
                             tau)
        out[win["bg"]] = 0.0
        out[win["d1"]], out[win["r2"]], out[win["d2"]], out[win["r3"]] = j1d, j1r, j2d, j2r
        for (a, b), loss in zip(shades, p[9:]):
            # a shaded stretch at the first temperature reaches only rise 2 (the long hold at the second temperature
            # lets it die away); one at the second temperature reaches decay 2 and rise 3
            if a < hold:
                d, ph, c, qq, after = d1, phi1, c0, q1, win["r2"]
            else:
                d, ph, c, qq, after = d2, phi2, c0b, q2, win["d2"] | win["r3"]
            sel = np.flatnonzero(after & (t > a))
            if sel.size:
                out[sel] -= loss * (step_flux(t[sel] - a, d, ph, c, qq, L, tau) - step_flux(t[sel] - b, d, ph, c, qq, L, tau))
        return pred + scale * out

    use = win["bg"] | win["d1"] | win["r2"] | win["d2"] | win["r3"]
    for sw in list(offs[:2]) + list(ons[1:3]):          # the readings within 5 s of a switch carry its electrical spike
        use &= ~((t >= sw) & (t < sw + 5.0))

    state = {"use": use.copy()}

    def resid(p):
        return (model(p) - cur)[state["use"]]

    smooth = median_filter(cur, size=15, mode="nearest")    # start values only: disturbances must not steer them
    b0 = np.median(cur[win["bg"]])
    j_ss1 = np.median(smooth[win["r2"]][-50:]) - b0
    rise = smooth[win["r2"]] - b0
    t_half = (t[win["r2"]] - ons[1])[np.argmax(rise >= 0.5 * j_ss1)]
    d1 = 0.138 * L * L / max(t_half, 1.0)
    j_ss2 = np.median(smooth[win["r3"]][-50:]) - b0
    rise3 = smooth[win["r3"]] - b0
    t_half3 = (t[win["r3"]] - ons[2])[np.argmax(rise3 >= 0.5 * j_ss2)]
    d2 = 0.138 * L * L / max(t_half3, 1.0)
    c0 = j_ss1 / scale * L / d1
    global N_TERMS                                      # 60 modes: ample after the 5 s left out at each switch
    N_TERMS = 60
    best = None
    q_top = 3 if exit_fit else 1e-12
    sigma = np.std(np.diff(cur[win["bg"]])) / math.sqrt(2) if win["bg"].sum() > 10 else 0.004
    n_sh = len(shades)
    tight = dict(ftol=1e-10, xtol=1e-10, gtol=1e-10)

    def refine(x, key, losses=True, **kw):
        """least_squares over the parameters that are free for this exit model (the others stay where x has them)."""
        lo, hi = bounds[key]
        free = np.array([True, True, True, c0_split, key == "finite", key == "finite", True, True, free_traps]
                        + [losses] * n_sh)
        idx = np.flatnonzero(free)
        x = np.array(x, float)

        def f(z):
            pp = x.copy()
            pp[idx] = z
            return resid(pp)

        sol = least_squares(f, np.clip(x[idx], lo[idx], hi[idx]), bounds=(lo[idx], hi[idx]), x_scale="jac", **kw)
        x[idx] = sol.x
        return SimpleNamespace(x=x, cost=sol.cost)

    trap0 = math.log(max(n_sheet, 1.0))
    bounds = {k: (np.array([-40, -40, -20, -20, 0, 0, -10, -10, math.log(0.05) if free_traps else trap0 - 1]
                           + [0.0] * n_sh, float),
                  np.array([-10, -10, 10, 10, q_top if k == "finite" else 1.0, q_top if k == "finite" else 1.0, 10, 10,
                            math.log(500.0) if free_traps else trap0 + 1] + [0.95] * n_sh, float))
              for k in ("sound", "finite")}
    # first without the shaded stretches: their readings and the after-effect (four breakthrough lags) are left out
    # and their losses held at zero; then every loss is freed and the whole record refitted from there
    calm = use.copy()
    for a, b in shades:
        lag = L * L / (6 * (d1 if a < hold else d2)) + tau_entry
        calm &= ~((t >= a) & (t < b + 4 * lag))
    state["use"] = calm
    fits = {}
    for key, q0 in [("sound", 0.0)] + ([("finite", 0.0), ("finite", 0.25)] if exit_fit else []):
        start = [math.log(d1 * phi1), math.log(d2 * phi2), math.log(max(c0 * phi1, 1e-6)),
                 math.log(max(c0 * phi1, 1e-6)), q0, q0, b0, 0.0, trap0] + [0.0] * n_sh
        fit = refine(start, key, losses=False, max_nfev=300, **tight)
        if key not in fits or fit.cost < fits[key].cost:
            fits[key] = fit
    state["use"] = use.copy()
    if n_sh:
        for key, fit in fits.items():
            fits[key] = refine(np.r_[fit.x[:9], [0.3] * n_sh], key, max_nfev=200, **tight)
    # disturbed stretches (sustained offsets, dips or bursts that the diffusion model cannot follow) must not pull the
    # fit: both exit models are refitted with a Cauchy loss (scale 2 noise units); a finite exit rate is kept only when
    # it lowers that robust misfit well beyond what the noise allows.  Residual windows of 15 readings whose mean or RMS
    # stand out of the noise are then left out (widened by a window on each side) and the plain fit is repeated, up to
    # three rounds, while the left-out share stays below 60 %.
    if not screen:                                      # shortcut: plain chi-square choice, nothing left out
        best = fits["sound"]
        if "finite" in fits and 2 * (fits["sound"].cost - fits["finite"].cost) / sigma ** 2 > 1000.0:
            best = fits["finite"]
    else:
        fs = 2 * sigma
        robust = {key: refine(fit.x, key, loss="cauchy", f_scale=fs, max_nfev=300) for key, fit in fits.items()}
        best_key = "sound"
        if "finite" in robust and 2 * (robust["sound"].cost - robust["finite"].cost) / fs ** 2 > 250.0:
            best_key = "finite"
        w = 15
        box = np.ones(w) / w
        best, flagged = robust[best_key], np.zeros(t.size, bool)
        for _ in range(3):
            r = np.where(use, model(best.x) - cur, 0.0)
            n_in = np.convolve(use.astype(float), box, "same")
            mean = np.convolve(r, box, "same") / np.maximum(n_in, 1e-9)
            rms = np.sqrt(np.convolve(r * r, box, "same") / np.maximum(n_in, 1e-9))
            bad = ((np.abs(mean) > 6 * sigma / math.sqrt(w)) | (rms > 4 * sigma)) & (n_in > 0.5)
            bad = np.convolve(bad.astype(float), np.ones(2 * w + 1), "same") > 0
            if not (bad & ~flagged).any() or (flagged | bad)[use].mean() > 0.6:
                break
            flagged |= bad
            state["use"] = use & ~flagged
            best = refine(best.x, best_key, max_nfev=300, **tight)
    ld1, ld2, lc0 = best.x[:3]
    d1, d2 = math.exp(ld1), math.exp(ld2)
    ea = R_GAS * math.log(d2 / d1) / (1 / T1 - 1 / T2) if ea_fixed is None else ea_fixed
    d25 = d1 * math.exp(-ea / R_GAS * (1 / 298.15 - 1 / T1))
    return {"diffusivity_25c_um2_per_s": d25 * 1e12, "subsurface_hydrogen_mol_per_m3": math.exp(lc0)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_membrane(sys.argv[1])))
