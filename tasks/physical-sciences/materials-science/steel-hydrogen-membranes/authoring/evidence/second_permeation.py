"""Second, separately written reduction of a permeation test (same call as the reference).

Instead of eigenfunction series it discretises the membrane by the method of lines (60 cells, finite-volume exit cell
with oxidation rate k, the sub-surface concentration as an extra state relaxing with the entry time constant, the
diffusion terms divided by the reversible-trap capacity) and advances each transient with a matrix exponential of the
constant-coefficient system, sampled on the data times by an eigen-decomposition of that matrix. The entry time
constant comes from one joint fit of the whole entry-potential record to the charging schedule. Parameters: D and
k L / D at each of the two temperatures, C0, a linear background and (where the sheet gives a binding energy) the
reversible-trap density, whose capacity 1 + (N_T / N_L) exp(E_b / R T) is tied to the two temperatures through the
sheet's E_b (the sheet's N_T only starts the fit), all fitted by Levenberg-Marquardt in log space; a
finite exit rate is kept only when it lowers chi-square by far more than the noise allows; D at 25 C from the
Arrhenius line through the two temperatures.  Bubbles shading the charged face show in the residual of that
potential fit as stretches lying beyond the charging value (21-point running median past 1 mV for 30 s or more); the
entry lag is refitted without them and each becomes a stretch of lowered entry target, its loss fitted, inside the
piecewise-constant schedule that the matrix exponentials follow (one in the hold changes the state decay 2 starts
from).  Disturbances lasting minutes are screened out between fits: readings
whose 21-point running median residual leaves the noise band (3 sigma), or whose own residual exceeds 8 sigma, are
dropped (widened by 21 readings) and the fits repeated until the kept set stops changing; start values come from a
31-point running median of the current, and the exit model is chosen on Cauchy-loss refits (scale 2 sigma) of the
plain fits.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import medfilt

FARADAY, GAS = 96485.33212, 8.314462618
CELLS = 60
N_LAT = 8.5e5


def system(d, bi, L, phi=1.0, tau=0.0):
    """dX/dt = A X + s * C0_target for cell averages X (plus the sub-surface concentration as the first state when
    tau > 0); exit flux g * C_last. Reversible traps divide the diffusion terms by their capacity phi."""
    h = L / CELLS
    k = bi * d / L
    o = 1 if tau > 0 else 0
    n = CELLS + o
    a = np.zeros((n, n))
    for i in range(CELLS):
        j = i + o
        if i > 0:
            a[j, j - 1] += d / h ** 2 / phi
            a[j, j] -= d / h ** 2 / phi
        if i < CELLS - 1:
            a[j, j + 1] += d / h ** 2 / phi
            a[j, j] -= d / h ** 2 / phi
    a[o, o] -= 2 * d / h ** 2 / phi                    # entry face at half a cell
    a[-1, -1] -= 1 / (h / (2 * d) + 1 / k) / h / phi   # exit: half-cell diffusion in series with oxidation
    s = np.zeros(n)
    if tau > 0:
        a[0, 0] = -1 / tau                              # sub-surface concentration relaxing towards the target
        a[o, 0] = 2 * d / h ** 2 / phi
        s[0] = 1 / tau
    else:
        s[0] = 2 * d / h ** 2 / phi
    g = 1 / (h / (2 * d) + 1 / k)                       # exit flux = g * C_last
    return a, s, g


_EIG = {}


def transient(a, s, g, c_start, entry, times):
    """Exit flux at `times` and the state at the last time, starting from c_start with the target held at `entry`."""
    key = a.tobytes()                       # the same system recurs within one prediction: decompose it once
    if key not in _EIG:
        if len(_EIG) > 16:
            _EIG.clear()
        _EIG[key] = np.linalg.eig(a)
    vals, vecs = _EIG[key]
    steady = -np.linalg.solve(a, s * entry)
    coef = np.linalg.solve(vecs, c_start - steady)
    flux = g * ((np.exp(np.outer(times, vals)) @ (vecs[-1] * coef)).real + steady[-1])
    last = (vecs @ (coef * np.exp(vals * times[-1]))).real + steady if times.size else c_start
    return flux, last


def schedule_response(t, switches, states, tau):
    """First-order response x(t) to the charging schedule (target 1 while charging, 0 otherwise)."""
    x = np.zeros(t.size)
    level = 0.0
    for k, (sw, st) in enumerate(zip(switches, states)):
        nxt = switches[k + 1] if k + 1 < len(switches) else np.inf
        sel = (t >= sw) & (t < nxt)
        x[sel] = st + (level - st) * np.exp(-(t[sel] - sw) / tau)
        level = st + (level - st) * math.exp(-min(nxt - sw, 1e7) / tau)
    return x


def entry_lag(t, pot, switches, states, ok=None):
    """One joint fit of the potential record (readings in `ok`), E = E_off + dE x(t), over the entry time constant.
    Returns the time constant and the residual of the whole record with its dE."""
    ok = np.ones(t.size, bool) if ok is None else ok

    def solve(lt):
        a = np.column_stack([np.ones(t.size), schedule_response(t, switches, states, math.exp(lt))])
        coef, *_ = np.linalg.lstsq(a[ok], pot[ok], rcond=None)
        return a @ coef - pot, coef[1]

    cost = lambda lt: float(np.sum(solve(lt)[0][ok] ** 2))
    grid = np.linspace(math.log(1e-2), math.log(3000.0), 60)
    k = int(np.argmin([cost(g) for g in grid]))
    lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]
    lt = minimize_scalar(cost, bounds=(lo, hi), method="bounded").x
    r, de = solve(lt)
    return math.exp(lt), -r, de


def shaded(t, resid, de, charging):
    """Stretches of charging where the potential lies beyond the fitted charging value: 21-point running median of the
    residual past 1 mV in the direction charging moves the potential, for 15 readings (30 s) or more."""
    beyond = (medfilt(resid, 21) * np.sign(de) > 1.0) & charging
    edges = np.flatnonzero(np.diff(np.r_[0, beyond.astype(int), 0]))
    step = t[1] - t[0]
    return [(t[i0] - step / 2, t[i1 - 1] + step / 2) for i0, i1 in zip(edges[::2], edges[1::2]) if i1 - i0 >= 15]


def run_rise(a, s, g, c_start, level, t0, times, cuts):
    """Exit flux at `times` (from t0) and the end state, the target at `level` except (start, end, loss) stretches."""
    flux = np.zeros(times.size)
    marks = [t0] + sorted({x for a0, b0, _ in cuts for x in (a0, b0) if x > t0}) + [np.inf]
    state = c_start
    for lo, hi in zip(marks[:-1], marks[1:]):
        loss = sum(l for a0, b0, l in cuts if a0 <= lo < b0)
        sel = (times >= lo - t0) & (times < hi - t0)
        if sel.any():
            flux[sel], _ = transient(a, s, g, state, level * (1 - loss), times[sel] - (lo - t0))
        if np.isfinite(hi):
            _, state = transient(a, s, g, state, level * (1 - loss), np.array([hi - lo]))
    return flux


def reduce_membrane(folder):
    folder = Path(folder)
    with open(folder / "transient.csv", newline="") as fh:
        data = np.array([[float(v) for v in r.values()] for r in csv.DictReader(fh)])
    t, cur, temp, pot = data.T
    with open(folder / "events.csv", newline="") as fh:
        ev = [(float(r["time_s"]), r["event"]) for r in csv.DictReader(fh)]
    on = [x for x, e in ev if e.startswith("charging on")]
    off = [x for x, e in ev if e.startswith("charging off")]
    change = [x for x, e in ev if e.startswith("thermostat")][0]
    thick = np.mean(np.loadtxt(folder / "micrometer.csv", delimiter=",", skiprows=1)[:, 1]) * 1e-3
    sheet = json.loads((folder / "specimen.json").read_text())
    area = sheet["exposed_area_cm2"] * 1e-4
    nt = sheet.get("reversible_trap_density_mol_per_m3", 0.0)
    eb = sheet.get("trap_binding_energy_kj_per_mol", 0.0) * 1e3
    parts = [(t >= off[0]) & (t < on[1]), (t >= on[1]) & (t < change), (t >= off[1]) & (t < on[2]), t >= on[2]]
    base = t < on[0]
    temps = [temp[parts[0] | parts[1]].mean() + 273.15, temp[parts[2] | parts[3]].mean() + 273.15]
    phis = [1 + nt / N_LAT * math.exp(eb / (GAS * tk)) for tk in temps]
    sw = sorted([(x, 1.0) for x in on] + [(x, 0.0) for x in off])
    tau, resid, de = entry_lag(t, pot, [x for x, _ in sw], [s for _, s in sw])
    charging = ((t >= on[1]) & (t < off[1])) | (t >= on[2])
    for x in on[1:3]:
        charging &= t > x + 5 * tau + 60.0
    cuts = shaded(t, resid, de, charging)
    if cuts:
        ok = np.ones(t.size, bool)
        for a0, b0 in cuts:
            ok &= ~((t >= a0) & (t < b0))
        tau = entry_lag(t, pot, [x for x, _ in sw], [s for _, s in sw], ok)[0]
    tau = tau if tau > 0.5 else 0.0
    gain = FARADAY * area * 1e6
    stride = max(1, t.size // 2000)
    keep = np.arange(t.size) % stride == 0
    for x in off[:2] + on[1:3]:
        keep &= ~((t >= x) & (t < x + 5.0))

    def predict(p, bi_fixed=None):
        out = p[5] + p[6] * t / 3600
        flux = np.zeros(t.size)
        for level, (i_dec, i_rise, t_off, t_on) in enumerate([(0, 1, off[0], on[1]), (2, 3, off[1], on[2])]):
            d = math.exp(p[level])
            bi = bi_fixed if bi_fixed is not None else math.exp(p[2 + level])
            phi = 1 + math.exp(p[7]) / N_LAT * math.exp(eb / (GAS * temps[level])) if eb > 0 else phis[level]
            a, s, g = system(d, bi, thick, phi, tau)
            c0 = math.exp(p[4])
            c_on = -np.linalg.solve(a, s * c0)
            mine = [(a0, b0, p[8 + k]) for k, (a0, b0) in enumerate(cuts)
                    if (a0 < change) == (level == 0)]
            held = sorted(c for c in mine if c[0] < t_off)              # in the hold: the decay starts elsewhere
            if held:
                state, t_now = c_on, held[0][0]
                for k, (a0, b0, l) in enumerate(held):
                    nxt = held[k + 1][0] if k + 1 < len(held) else t_off
                    _, state = transient(a, s, g, state, c0 * (1 - l), np.array([b0 - a0]))
                    _, state = transient(a, s, g, state, c0, np.array([nxt - b0]))
                c_on = state
            sel_d, sel_r = parts[i_dec] & keep, parts[i_rise] & keep
            j_dec, _ = transient(a, s, g, c_on, 0.0, t[sel_d] - t_off)
            _, c_switch = transient(a, s, g, c_on, 0.0, np.array([t_on - t_off]))
            j_rise = run_rise(a, s, g, c_switch, c0, t_on, t[sel_r] - t_on, [c for c in mine if c[0] >= t_on])
            flux[sel_d], flux[sel_r] = j_dec, j_rise
        return out + gain * flux

    mask = (base | parts[0] | parts[1] | parts[2] | parts[3]) & keep
    calm = medfilt(cur, 31)                              # start values from a running median, not from raw readings
    b0 = np.median(cur[base])
    jss = [np.median(calm[parts[1]][-40:]) - b0, np.median(calm[parts[3]][-40:]) - b0]
    guess_d = []
    for idx, t_on in ((1, on[1]), (3, on[2])):
        y = calm[parts[idx]] - b0
        tt = t[parts[idx]] - t_on
        half = tt[np.argmax(y >= 0.5 * np.median(y[-40:]))]
        guess_d.append(0.138 * thick ** 2 / max(half, 1.0))
    c0 = jss[0] / gain * thick / guess_d[0]
    sigma = np.std(np.diff(cur[base])) / math.sqrt(2)
    p0 = [math.log(guess_d[0] * phis[0]), math.log(guess_d[1] * phis[1]), 0.0, 0.0,
          math.log(max(c0 * phis[0], 1e-6)), b0, 0.0, math.log(max(nt, 1.0))] + [0.2] * len(cuts)
    idx_fit = np.flatnonzero(mask)

    def fit_on(use, start, finite):
        fun = (lambda p: (predict(p) - cur)[use]) if finite else (lambda p: (predict(p, 1e6) - cur)[use])
        return least_squares(fun, start, method="lm", x_scale="jac", max_nfev=200)

    def screen(p, finite):
        # readings whose 21-point running median residual leaves the noise (3 sigma) or whose own residual exceeds 8
        # sigma are dropped, widened by 21 readings: disturbances that last minutes do not belong to the diffusion
        r = (predict(p) if finite else predict(p, 1e6)) - cur
        rf = r[idx_fit]
        med = medfilt(rf, 21)
        bad = (np.abs(med) > 3 * sigma) | (np.abs(rf) > 8 * sigma)
        bad = np.convolve(bad.astype(float), np.ones(21), "same") > 0
        use = mask.copy()
        use[idx_fit[bad]] = False
        return use

    use = mask
    lo = np.full(8 + len(cuts), -np.inf)
    hi = np.full(8 + len(cuts), np.inf)
    lo[8:], hi[8:] = 0.0, 0.95
    lo[2:4], hi[2:4] = math.log(0.1), math.log(1e4)
    lo[0:2], hi[0:2] = -40.0, -10.0
    lo[4], hi[4] = -20.0, 10.0
    lo[7], hi[7] = -10.0, 10.0
    # plain fits of both exit models, then Cauchy-loss refits (scale 2 sigma) from them, so that disturbances do not
    # decide the exit model; the finite exit rate is kept when it lowers the robust misfit far beyond the noise
    sound = fit_on(use, p0, False)
    cands = [(False, sound.x)]
    for lbi in (math.log(50.0), math.log(4.0)):
        q = np.array(sound.x)
        q[2] = q[3] = lbi
        cands.append((True, fit_on(use, q, True).x))
    scored = []
    for fin, x in cands:
        fun = (lambda p: (predict(p) - cur)[mask]) if fin else (lambda p: (predict(p, 1e6) - cur)[mask])
        try:
            r = least_squares(fun, np.clip(x, lo + 1e-9, hi - 1e-9), method="trf", loss="cauchy", f_scale=2 * sigma,
                              x_scale="jac", bounds=(lo, hi), max_nfev=100)
        except (OverflowError, np.linalg.LinAlgError, ValueError):
            continue
        scored.append((fin, r))
    best_sound = scored[0][1]
    best, finite = best_sound, False
    for fin, r in scored[1:]:
        if 2 * (best_sound.cost - r.cost) / (2 * sigma) ** 2 > 60.0 and r.cost < best.cost:
            best, finite = r, True
    for _ in range(4):
        new = screen(best.x, finite)
        if new.sum() < 0.4 * mask.sum() or np.array_equal(new, use):
            break
        use = new
        best = fit_on(use, best.x, finite)
    d1, d2 = math.exp(best.x[0]), math.exp(best.x[1])
    ea = GAS * math.log(d2 / d1) / (1 / temps[0] - 1 / temps[1])
    d25 = d1 * math.exp(-ea / GAS * (1 / 298.15 - 1 / temps[0]))
    return {"diffusivity_25c_um2_per_s": d25 * 1e12, "subsurface_hydrogen_mol_per_m3": math.exp(best.x[4])}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_membrane(sys.argv[1])))
