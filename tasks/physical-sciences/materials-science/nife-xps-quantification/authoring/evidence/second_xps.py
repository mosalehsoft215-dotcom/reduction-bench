"""Second, independent reduction of one electrode (same notes, different numerics).

    python second_xps.py <electrode directory> <result.json>

- Parsing with NumPy only; regions split by label.
- Line shapes: every component's intrinsic width combined in quadrature with the analyser broadening 0.012 Ep and,
  for the charged film, with a region's uneven-charging broadening 2.3548 kappa |q| at the region's mean charge (the
  grounded graphite has none).
- Uneven charging: the extra broadening of the closing C 1s region, recorded once the charge has settled, is scanned
  on a 0.05 eV grid and then refined, the position refitted at each trial (graphite unbroadened, at the pass energy
  the graphite chose), and gives
  kappa = broadening / (2.3548 |q|) (zero when |q| < 0.5 eV, where it cannot matter); the C 1s areas and positions,
  the build-up, and the O 1s, Ni 2p3/2 and Fe 2p3/2 pass energies and fits are then done with the broadening of each
  region.
- Energy scale: in each C 1s region one least-squares fit of the main C 1s template and the grounded graphite with
  both positions free, from a grid of starts and at Ep = 20 and 50 eV (the lower misfit wins); the graphite position
  over 284.40 eV gives the scale gain, the mean of the two regions is used, and every axis is divided by (1 + gain).
- Pass energies actually used: for C 1s, Ni 2p3/2 and O 1s the one of 20 or 50 eV whose broadened template has the
  smaller variable-projection misfit over a 0.05 eV position scan refined by a bounded search; for Fe 2p3/2 the one whose full window fit (Auger
  included) has the smaller cost.  The logged pass energies and the logged scale gain are not used.
- Positions (C 1s, O 1s): variable projection, a 0.03 eV scan of the main position with heights, a straight
  background, the Shirley share and (on carbon paper) the grounded graphitic profile linear at each trial, then a
  bounded one-dimensional refinement.
- Charge: the region-averaged charges of the two C 1s and the lattice O 1s fix the first-order build-up, solved by a
  scan of its time constant on a log grid with the two levels linear.
- Each metal window: one nonlinear least-squares fit of position, height, constant, slope and Shirley factor
  together, the Shirley step built from the model peaks themselves (Auger included in the Fe window, its height fixed
  by the Ni result).
- Areas: the fitted template's area inside the recorded window; rates per second of recording.
- Overlayer: the adventitious C 1s (opening and closing) and O 1s area rates from linear fits at the found positions;
  the thickness at both ends, linear in time between them, from the two balances
  A_C = 0.85 (1 - exp(-d / lambda_C)) sum_X A_X exp(d(t_X) / lambda_X), one per C 1s region, solved together (fsolve)
  from the single-thickness fixed-point solution, when the two C 1s rates differ by more than 16 %.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import least_squares, minimize_scalar

HV, DE = 1486.6, 0.05
SENS = {"Ni2p3": 14.6, "Fe2p3": 10.8}             # C 1s 1.0 and O 1s 2.93 enter the overlayer step
XEXP = {50: 0.6, 20: 0.85}
AUG = (712.6, 3.35, 0.07)
LINES = {"C1s": ((0.0, 1.15, 1.0), (1.5, 1.25, 0.18), (3.8, 1.25, 0.08)),
         "Ni2p3": ((0.0, 2.0, 1.0), (1.8, 2.5, 0.45), (6.1, 3.95, 0.35))}
FE_SAT = ((713.2, 2.95, 0.35), (718.8, 3.45, 0.12))
FE_MAIN = 2.5
PE = [50.0]                                   # pass energy of the region being modelled (set per region)
EXTRA = [0.0]                                 # uneven-charging FWHM of the region being modelled (charged film only)
SWEEP = [((0.0, 0.0),)]                       # (charge offset from the region mean, uneven-charging FWHM) per instant


def profile(e, centre, fwhm):
    """Line of the charged film: averaged over the instants of the region's charge sweep, each with its own
    uneven-charging broadening, on top of the intrinsic width and the analyser broadening."""
    out = 0.0
    for dq, bb in SWEEP[0]:
        f = np.sqrt(fwhm * fwhm + (0.012 * PE[0]) ** 2 + EXTRA[0] ** 2 + bb * bb)
        z = (e - centre - dq) / f
        out = out + 0.7 * np.exp(-4 * np.log(2) * z * z) + 0.3 / (1 + 4 * z * z)
    return out / len(SWEEP[0])


def profile_g(e, centre, fwhm):
    """The grounded graphite: no charge, so neither sweep nor uneven-charging broadening."""
    keep, sweep, EXTRA[0], SWEEP[0] = EXTRA[0], SWEEP[0], 0.0, ((0.0, 0.0),)
    out = profile(e, centre, fwhm)
    EXTRA[0], SWEEP[0] = keep, sweep
    return out


def peaks(region, e, pos, shift):
    if region == "Fe2p3":
        y = profile(e, pos, FE_MAIN)
        for c, w, h in FE_SAT:
            y = y + h * profile(e, c + shift, w)
        return y
    return sum(h * profile(e, pos + off, w) for off, w, h in LINES[region])


def template_area(region, pos, shift):
    grid = np.arange(pos - 400, pos + 400, 0.01)
    return np.trapezoid(peaks(region, grid, pos, shift), grid)


def model(p, region, e, shift, aug_curve):
    pos, height, c0, c1, ks = p
    pk = height * peaks(region, e, pos, shift) + aug_curve
    return pk, c0 + c1 * (e - e[0]) + ks * np.cumsum(pk) * DE


def nonlinear_fit(region, e, y, start, shift=0.0, aug_curve=None):
    aug_curve = np.zeros_like(e) if aug_curve is None else aug_curve
    w = 1 / np.sqrt(np.maximum(y, 1.0))
    base = np.polyfit(e - e[0], np.r_[y[:20], y[-20:]].mean() + 0 * e, 0)[0]
    height0 = max(y.max() - base, 1.0)
    res = lambda p: (sum(model(p, region, e, shift, aug_curve)) - y) * w
    sol = least_squares(res, [start, height0, y[:10].mean(), 0.0, 0.0], x_scale=[0.1, height0, 100, 1, 0.01],
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
    nonlinear_fit.cost = sol.cost
    return sol.x


O1S = ((0.0, 1.05, 1.0), (1.7, 1.5, 0.42))
GRAPH = (284.40, 0.62)


def line_position(e, y, comps, extras=()):
    """Main-component position by variable projection: a 0.01 eV scan, then a bounded 1-D refinement; heights,
    a straight background, the Shirley share and any extra fixed-position profiles are linear at each trial."""
    def cost(pos):
        s = sum(h * profile(e, pos + off, w) for off, w, h in comps)
        cols = [s, np.cumsum(s) * DE, np.ones_like(e), e - e[0]]
        for x in extras:
            cols += [x, np.cumsum(x) * DE]
        a = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(a, y, rcond=None)
        return np.sum((a @ coef - y) ** 2 / np.maximum(y, 1))
    grid = np.arange(e[0] + 2, e[-1] - 4, 0.03)
    k = int(np.argmin([cost(x) for x in grid]))
    return minimize_scalar(cost, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]), method="bounded",
                           options={"xatol": 1e-6}).x


def window_area(e, y, comps, pos, extras=()):
    """Template area inside the window at a fixed main position (linear fit of height, background and steps)."""
    s = sum(h * profile(e, pos + off, w) for off, w, h in comps)
    cols = [s, np.cumsum(s) * DE, np.ones_like(e), e - e[0]]
    for x in extras:
        cols += [x, np.cumsum(x) * DE]
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)
    return coef[0] * np.sum(s) * DE


def scan_cost(e, y, comps, extras=()):
    """Smallest variable-projection misfit of a template over a 0.02 eV scan of its main position."""
    def cost(pos):
        sh = sum(h * profile(e, pos + off, w) for off, w, h in comps)
        cols = [sh, np.cumsum(sh) * DE, np.ones_like(e), e - e[0]]
        for x in extras:
            cols += [x, np.cumsum(x) * DE]
        a = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(a, y, rcond=None)
        return np.sum((a @ coef - y) ** 2 / np.maximum(y, 1))
    grid = np.arange(e[0] + 2, e[-1] - 4, 0.05)
    k = int(np.argmin([cost(x) for x in grid]))
    return minimize_scalar(cost, bounds=(grid[max(k - 1, 0)], grid[min(k + 1, grid.size - 1)]), method="bounded",
                           options={"xatol": 1e-3}).fun


def used_pass(e, y, comps, extras_fn=lambda e: ()):
    """The pass energy (20 or 50 eV) whose analyser broadening fits the region's line shapes better."""
    best = None
    for pe in (20.0, 50.0):
        PE[0] = pe
        c = scan_cost(e, y, comps, extras_fn(e))
        if best is None or c < best[0]:
            best = (c, pe)
    PE[0] = best[1]
    return best[1]


def graphite_gain(e, y):
    """Scale gain from the grounded graphite of carbon paper: free main C 1s and graphite positions in one fit, at
    each candidate pass energy; the better-fitting pass energy wins."""
    return min((graphite_fit(e, y, pe) for pe in (20.0, 50.0)), key=lambda r: r[1])[0]


def graphite_fit(e, y, pe):
    PE[0] = pe
    def res(p):
        main, gpos, hm, hg, c0, c1, ks = p
        pk = hm * sum(h * profile(e, main + off, w) for off, w, h in LINES["C1s"]) + hg * profile_g(e, gpos, GRAPH[1])
        return (pk + c0 + c1 * (e - e[0]) + ks * np.cumsum(pk) * DE - y) / np.sqrt(np.maximum(y, 1.0))
    top = y.max() - y[:10].mean()
    best = None
    for g0 in GRAPH[0] * (1 + np.arange(-0.003, 0.00301, 0.0005)):
        for m0 in e[np.argmax(y)] + np.array([-1.5, -0.5, 0.5, 1.5]):
            sol = least_squares(res, [m0, g0, top, top / 2, y[:10].mean(), 0.0, 0.0], max_nfev=200,
                                x_scale=[0.1, 0.1, top, top, 100, 1, 0.01])
            if best is None or sol.cost < best.cost:
                best = sol
    return best.x[1] / GRAPH[0] - 1.0, best.cost


def att(ke):
    return 0.0175 * ke ** 0.75


def unattenuate(a_c, a_ox, ke_c, ke_ox):
    """Overlayer thickness by fixed-point iteration on d = -lam_C ln(1 - A_C / (0.85 k(d))), k(d) = sum A_X e^(d/lam_X)."""
    d = 0.0
    for _ in range(200):
        k = sum(a * np.exp(d / att(ke)) for a, ke in zip(a_ox, ke_ox))
        ratio = min(a_c / (0.85 * k), 0.999)
        d_new = -att(ke_c) * np.log(1 - ratio)
        if abs(d_new - d) < 1e-10:
            break
        d = 0.5 * d + 0.5 * d_new
    return d


def build_up(obs, spans):
    """q(t) = qi + (q0 - qi) exp(-t / tau): for each tau on a log grid q0, qi are linear; the best tau wins."""
    best = None
    for tau in np.geomspace(30.0, 3e4, 400):
        rows = []
        for (t0, t1) in spans:
            f = tau * (np.exp(-t0 / tau) - np.exp(-t1 / tau)) / (t1 - t0)
            rows.append([f, 1 - f])
        a = np.array(rows)
        coef, *_ = np.linalg.lstsq(a, np.array(obs), rcond=None)
        r = np.sum((a @ coef - obs) ** 2)
        if best is None or r < best[0] - 1e-12:
            best = (r, tau, coef)
    _, tau, (q0, qi) = best
    return lambda t0, t1: q0 * tau * (np.exp(-t0 / tau) - np.exp(-t1 / tau)) / (t1 - t0) + qi * (
        1 - tau * (np.exp(-t0 / tau) - np.exp(-t1 / tau)) / (t1 - t0))


def integrated(region, e, y, p, shift, aug_curve):
    # the template's area inside the recorded window, as the laboratory's areas are defined
    return p[1] * np.sum(peaks(region, e, p[0], shift)) * DE


def trans(pe, ke):
    return pe / 50 * (ke / 1000) ** (-XEXP[pe])


def run(folder):
    folder = Path(folder)
    sched = json.loads((folder / "schedule.json").read_text())
    raw = np.genfromtxt(folder / "regions.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    reg = {}
    for lab in set(raw["region"].tolist()):
        m = raw["region"] == lab
        reg[lab] = (raw["binding_energy_ev"][m].astype(float), raw["counts"][m].astype(float))
    info = {r["region"]: r for r in sched["regions"]}
    tmid = {k: (v["start_s"] + v["end_s"]) / 2 for k, v in info.items()}
    span = {k: v["end_s"] - v["start_s"] for k, v in info.items()}
    paper = sched.get("substrate") == "carbon paper"
    # energy scale: the grounded graphite fixes the scale gain; every axis is put back on (E + q)
    if paper:
        gains = []
        for lab in ("C1s", "C1s_end"):
            e, y = reg[lab]
            gains.append(graphite_gain(e, y))
        gain = float(np.mean(gains))
    else:
        gain = sched["scale_gain_ppm"] * 1e-6
    reg = {k: (e / (1 + gain), y) for k, (e, y) in reg.items()}
    fe_shape = ((0.0, FE_MAIN, 1.0), (2.3, FE_SAT[0][1], 0.35), (7.9, FE_SAT[1][1], 0.12))
    shapes = {"C1s": LINES["C1s"], "C1s_end": LINES["C1s"], "Ni2p3": LINES["Ni2p3"], "Fe2p3": fe_shape, "O1s": O1S}
    pe_used = {}
    for lab in ("C1s", "C1s_end", "O1s"):
        e, y = reg[lab]
        extra = (lambda x: (profile_g(x, *GRAPH),)) if (paper and lab.startswith("C1s")) else (lambda x: ())
        pe_used[lab] = used_pass(e, y, shapes[lab], extra)
    ref_span = [(info[lab]["start_s"], info[lab]["end_s"]) for lab in ("C1s", "C1s_end", "O1s")]
    kappa, q_mean = 0.0, None

    def set_sweep(k):
        """The charge sweep of region k at six instants (one when settled), with kappa's broadening at each."""
        if q_mean is None:
            SWEEP[0] = ((0.0, 0.0),)
            return
        t0, t1 = info[k]["start_s"], info[k]["end_s"]
        ts = t0 + (np.arange(6) + 0.5) * (t1 - t0) / 6
        qs = np.array([q_mean(x - 0.05, x + 0.05) for x in ts])
        if np.ptp(qs) <= 0.02:          # a settled charge: one instant is enough
            qs = qs[:1] * 0 + qs.mean()
        SWEEP[0] = tuple((float(q - qs.mean()), float(2.3548 * kappa * abs(q))) for q in qs)

    for rnd in range(2):
        ref_q, c_rates, c_broad = [], {}, {}
        for lab in ("C1s", "C1s_end"):
            e, y = reg[lab]
            PE[0] = pe_used[lab]
            EXTRA[0] = 0.0
            set_sweep(lab)
            extra = (profile_g(e, *GRAPH),) if paper else ()
            pos = line_position(e, y, LINES["C1s"], extra)
            c_rates[lab] = window_area(e, y, LINES["C1s"], pos, extra) / span[lab] / trans(pe_used[lab], HV - 284.8)
            ref_q.append(pos - 284.8)
            if rnd == 0 and lab == "C1s_end":    # its extra broadening, scanned at the found position
                def at_pos(bb, x0):
                    EXTRA[0] = bb
                    sh = sum(h * profile(e, x0 + off, w) for off, w, h in LINES["C1s"])
                    cols = [sh, np.cumsum(sh) * DE, np.ones_like(e), e - e[0]]
                    for x in extra:
                        cols += [x, np.cumsum(x) * DE]
                    a_ = np.column_stack(cols)
                    coef, *_ = np.linalg.lstsq(a_, y, rcond=None)
                    return np.sum((a_ @ coef - y) ** 2 / np.maximum(y, 1))

                def misfit(bb):         # the position refitted at each trial broadening
                    return minimize_scalar(lambda x0: at_pos(bb, x0), bounds=(pos - 0.4, pos + 0.4),
                                           method="bounded", options={"xatol": 1e-4}).fun
                coarse = np.arange(0.0, 1.61, 0.05)
                k = int(np.argmin([misfit(bb) for bb in coarse]))
                c_broad[lab] = float(minimize_scalar(misfit, bounds=(coarse[max(k - 1, 0)], coarse[min(k + 1, coarse.size - 1)]),
                                                     method="bounded", options={"xatol": 1e-3}).x)
        e, y = reg["O1s"]
        EXTRA[0] = 0.0
        set_sweep("O1s")
        if q_mean is not None:
            pe_used["O1s"] = used_pass(e, y, O1S)
        PE[0] = pe_used["O1s"]
        o_pos = line_position(e, y, O1S)
        o_rate = window_area(e, y, O1S, o_pos) / span["O1s"] / trans(pe_used["O1s"], HV - 529.70)
        ref_q.append(o_pos - 529.70)
        q_mean = build_up(np.array(ref_q), ref_span)
        if rnd == 0:                    # the closing C 1s, recorded once the charge has settled, sets kappa
            q_end = abs(q_mean(*ref_span[1]))
            kappa = c_broad["C1s_end"] / (2.3548 * q_end) if q_end > 0.5 else 0.0
    shift = lambda k: q_mean(info[k]["start_s"], info[k]["end_s"])

    e, y = reg["Ni2p3"]
    set_sweep("Ni2p3")
    pe_used["Ni2p3"] = used_pass(e, y, shapes["Ni2p3"])
    PE[0] = pe_used["Ni2p3"]
    p = nonlinear_fit("Ni2p3", e, y, e[np.argmax(y)])
    ni_true = p[0] - shift("Ni2p3")
    ni_area = integrated("Ni2p3", e, y, p, 0.0, np.zeros_like(e))
    ni_rate0 = ni_area / span["Ni2p3"] / trans(pe_used["Ni2p3"], HV - ni_true)

    e, y = reg["Fe2p3"]
    sf = shift("Fe2p3")
    set_sweep("Fe2p3")
    full = np.trapezoid(profile(np.arange(AUG[0] - 400, AUG[0] + 400, 0.01), AUG[0], AUG[1]), dx=0.01)
    trials = []
    for pe_fe in (20.0, 50.0):          # the Fe window's pass energy: the full Fe model (Auger included) at each
        PE[0] = pe_fe
        full = np.trapezoid(profile(np.arange(AUG[0] - 400, AUG[0] + 400, 0.01), AUG[0], AUG[1]), dx=0.01)
        a_area = AUG[2] * ni_rate0 * trans(pe_fe, HV - AUG[0]) * span["Fe2p3"]
        aug_curve = profile(e, AUG[0] + sf, AUG[1]) * a_area / full
        p = nonlinear_fit("Fe2p3", e, y, 710.9 + sf, sf, aug_curve)
        trials.append((nonlinear_fit.cost, pe_fe, p, aug_curve))
    _, pe_fe, p, aug_curve = min(trials, key=lambda r: r[0])
    PE[0] = pe_fe
    fe_area = integrated("Fe2p3", e, y, p, sf, aug_curve)
    fe_rate0 = fe_area / span["Fe2p3"] / trans(pe_fe, HV - (p[0] - sf))
    n_fe, n_ni = fe_rate0 / SENS["Fe2p3"], ni_rate0 / SENS["Ni2p3"]
    ke_ni, ke_fe = HV - ni_true, HV - (p[0] - sf)
    # the overlayer may thicken during the session: a thickness at the opening and one at the closing C 1s, linear in
    # time between them, from the two carbon balances solved together (Powell hybrid, scipy fsolve)
    t_a, t_b = tmid["C1s"], tmid["C1s_end"]
    ox = [(n_ni, ke_ni, tmid["Ni2p3"]), (n_fe, ke_fe, tmid["Fe2p3"]), (o_rate / 2.93, HV - 529.70, tmid["O1s"])]
    da = db = unattenuate(c_rates["C1s"], [a for a, _, _ in ox], HV - 284.8, [k for _, k, _ in ox])
    at = lambda tt: da + (db - da) * (tt - t_a) / (t_b - t_a)
    lam_c = att(HV - 284.8)

    def gap(dd):
        f = lambda tt: dd[0] + (dd[1] - dd[0]) * (tt - t_a) / (t_b - t_a)
        tot = sum(a * np.exp(f(tt) / att(k)) for a, k, tt in ox)
        return [0.85 * (1 - np.exp(-dd[0] / lam_c)) * tot - c_rates["C1s"],
                0.85 * (1 - np.exp(-dd[1] / lam_c)) * tot - c_rates["C1s_end"]]
    if abs(c_rates["C1s_end"] / c_rates["C1s"] - 1) > 0.16:      # the carbon grew: solve both ends
        da, db = fsolve(gap, [da, db], xtol=1e-12)
    n_ni, n_fe = n_ni * np.exp(at(tmid["Ni2p3"]) / att(ke_ni)), n_fe * np.exp(at(tmid["Fe2p3"]) / att(ke_fe))
    return {"ni2p3_binding_ev": float(ni_true), "fe_fraction": float(n_fe / (n_fe + n_ni))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
