"""Reduces the XPS narrow regions of one electrode to the charge-corrected Ni 2p3/2 binding energy and the
Fe / (Ni + Fe) atomic fraction.

    python3 xps_reduce.py <electrode directory> <result.json>

1. Templates: every region is fitted with its component template from the notes (GL(30) shapes at fixed offsets
   and height ratios, one free position and one free scale) on a constant, a linear slope and a Shirley-type step
   proportional to the peak intensity accumulated from the low binding-energy end.  Each component's FWHM is its
   intrinsic width combined in quadrature with the analyser broadening 0.012 Ep and, for the charged film (not the
   grounded graphite), with the uneven-charging broadening 2.3548 kappa |q| at the instant's charge.
2. Pass energy: each region is fitted at both of the analyser's pass energies (20 and 50 eV); the one with the lower
   chi-square is the pass energy the region was recorded at, and it sets the widths, the transmission and the Auger
   intensity of that region.
3. Scale gain: the paper's graphitic carbon is grounded, so its fitted position in the two C 1s regions (position
   free, alongside the adventitious template) is (1 + g) 284.40 eV; g converts the axis to E + q = axis / (1 + g).
4. Charge: on that scale the adventitious C 1s (284.8 eV) of the first and last regions and the lattice O 1s
   (529.70 eV) give three region-averaged charges, which fix the build-up q(t) = q_inf + (q_0 - q_inf) exp(-t/tau);
   the Ni and Fe regions use its average over their own clock times.
4b. Uneven charging: with the build-up known, kappa is the value that best fits both C 1s regions (pass energy from
   the graphite, which does not charge, the adventitious template broadened by 2.3548 kappa |q(t)|); the scale, the
   build-up and kappa are found again with that broadening (two rounds), every later fit carries it, and the O 1s
   pass energy is chosen again with it.
5. Ni 2p3/2 binding energy: the fitted main position on the corrected scale less the Ni region's average charge.
6. Composition: the Ni LMM Auger line enters the Fe window with 0.07 x (Ni 2p3/2 area rate / T_Ni) x T_Auger at the
   Fe region's pass energy; template area rates over sensitivity factor and transmission give each element's
   attenuated amount; the overlayer thickness d solves 0.85 (1 - exp(-d/lambda_C)) sum_X A_X exp(d/lambda_X) = A_C,
   and Ni and Fe are multiplied by exp(d/lambda) before Fe / (Ni + Fe).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

HV, STEP = 1486.6, 0.05
RSF = {"C1s": 1.0, "Ni2p3": 14.6, "Fe2p3": 10.8, "O1s": 2.93}
T_EXP = {50: 0.6, 20: 0.85}
RES = 0.012                                  # analyser broadening, eV of FWHM per eV of pass energy
AUGER_BE, AUGER_W, AUGER_RATIO = 712.6, 3.35, 0.07
GRAPHITE_BE, GRAPHITE_W = 284.40, 0.62
O_LATTICE, RHO_C = 529.70, 0.85
# template: (offset from the main component eV, intrinsic FWHM eV, height relative to the main component)
TEMPLATE = {"C1s": [(0.0, 1.15, 1.0), (1.5, 1.25, 0.18), (3.8, 1.25, 0.08)],
            "Ni2p3": [(0.0, 2.0, 1.0), (1.8, 2.5, 0.45), (6.1, 3.95, 0.35)],
            "Fe2p3": [(0.0, 2.5, 1.0), (None, 2.95, 0.35), (None, 3.45, 0.12)],
            "O1s": [(0.0, 1.05, 1.0), (1.7, 1.5, 0.42)]}
FIXED = {"Fe2p3": (713.2, 718.8)}            # true positions of the Fe satellites; they move with the charge


def gl(e, c, w, h=1.0):
    s = w / 2.3548
    return h * (0.7 * np.exp(-0.5 * ((e - c) / s) ** 2) + 0.3 / (1 + ((e - c) / (w / 2)) ** 2))


def width(w, pe):
    return float(np.hypot(w, RES * pe))


def transmission(pe, ke):
    return (pe / 50) * (ke / 1000) ** (-T_EXP[pe])


def shape(key, e, pos, charge, pe, smear=(0.0,), kq=None):
    """The region's template, averaged over the charge offsets `smear` the region swept through; kq gives the
    uneven-charging FWHM at each of those instants."""
    kq = (0.0,) * len(smear) if kq is None else kq
    out = np.zeros_like(e)
    for k, (off, w, h) in enumerate(TEMPLATE[key]):
        c = pos + off if off is not None else FIXED[key][k - 1] + charge
        out += sum(gl(e, c + dq, float(np.hypot(width(w, pe), b)), h) for dq, b in zip(smear, kq)) / len(smear)
    return out


def fit_region(key, e, y, step, pe, charge=0.0, extra=None, lin=(), free=None, z0=0.0, smear=(0.0,), hint=None,
               kq=None):
    """Template fit on the corrected scale e.  `extra`: a fixed additional profile (counts) with its own share of
    the step; `lin`: fixed-position profiles whose heights (and step shares) are linear; `free`: a callable
    z -> profile whose position z is a second nonlinear parameter.  Returns main position, template area,
    chi-square and z."""
    lo = e[0]
    target = y - (0 if extra is None else extra)
    weight = 1 / np.sqrt(np.clip(y, 1, None))

    def design(p):
        s = shape(key, e, p[0], charge, pe, smear, kq)
        cols = [s, np.cumsum(s) * step, np.ones_like(e), e - lo]
        if extra is not None:
            cols.append(np.cumsum(extra) * step)
        for f in list(lin) + ([free(p[1])] if free is not None else []):
            cols += [f, np.cumsum(f) * step]
        return np.column_stack(cols), s

    def resid(p):
        a, _ = design(p)
        coef = np.linalg.lstsq(a * weight[:, None], target * weight, rcond=None)[0]
        return (a @ coef - target) * weight

    trial = np.arange(e[0] + 1.0, e[-1] - 1.0, 0.2)
    if free is None and hint is not None:
        start = [hint]
    elif free is None:
        start = [trial[int(np.argmin([np.sum(resid([x]) ** 2) for x in trial]))]]
    else:                           # main position scanned at z0, then the free position within 0.6 eV of z0
        grid = [(x, z) for x in trial for z in z0 + np.arange(-0.6, 0.61, 0.1)] if hint is None else [hint]
        start = list(grid[int(np.argmin([np.sum(resid(p) ** 2) for p in grid]))])
    best = least_squares(resid, start, diff_step=1e-6)
    a, s = design(best.x)
    coef = np.linalg.lstsq(a * weight[:, None], target * weight, rcond=None)[0]
    return best.x[0], coef[0] * s.sum() * step, 2 * best.cost, (best.x[1] if free is not None else None)


def mean_charge(q0, qi, tau, t0, t1):
    """Average of q(t) = qi + (q0 - qi) exp(-t / tau) over [t0, t1]."""
    return qi + (q0 - qi) * tau * (np.exp(-t0 / tau) - np.exp(-t1 / tau)) / (t1 - t0)


def imfp(ke):
    """Attenuation length of the notes, nm, for a kinetic energy in eV."""
    return 0.0175 * ke ** 0.75


def overlayer_thickness(a_c, a_ox, ke_c, ke_ox):
    """Solve A_C = RHO_C (1 - exp(-d / lam_C)) sum_X A_X exp(d / lam_X) for d (nm)."""
    def f(d):
        return RHO_C * (1 - np.exp(-d / imfp(ke_c))) * sum(a * np.exp(d / imfp(k)) for a, k in zip(a_ox, ke_ox)) - a_c
    return 0.0 if f(0.0) >= 0 else brentq(f, 0.0, 20.0)


def reduce_electrode(folder, gain="graphite", passes="widths", drift="build-up", own_pass=True, auger="own",
                     rates=True, substrate=True, overlayer=True, growth=True, uneven=True):
    folder = Path(folder)
    head = json.loads((folder / "schedule.json").read_text())
    data = {}
    with open(folder / "regions.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            data.setdefault(r["region"], []).append((float(r["binding_energy_ev"]), float(r["counts"])))
    meta = {m["region"]: m for m in head["regions"]}
    arr = {k: np.array(v) for k, v in data.items()}
    span = {k: (m["start_s"], m["end_s"]) for k, m in meta.items()}
    dur = {k: (m["end_s"] - m["start_s"]) if rates else 520.0 for k, m in meta.items()}
    g_log = 1 + head["scale_gain_ppm"] * 1e-6

    chosen = {}

    def pick(label, fits):
        """fits: pass energy -> fit result; the plan's value, or the lower chi-square of the two (kept once found)."""
        if passes == "plan":
            chosen[label] = meta[label]["pass_energy_ev"]
        if label in chosen:
            return chosen[label], fits(chosen[label])
        got = {p: fits(p) for p in (20, 50)}
        chosen[label] = min(got, key=lambda k: got[k][2])
        return chosen[label], got[chosen[label]]

    def c1s(label, scale, sm, kq=None):
        e, y = arr[label][:, 0] / (scale or 1.0), arr[label][:, 1]
        st = STEP / (scale or 1.0)
        if not substrate:
            return pick(label, lambda p: fit_region("C1s", e, y, st, p, smear=sm, kq=kq))
        if scale is None:               # raw axis, graphite position free
            got = pick(label, lambda p: fit_region("C1s", e, y, STEP, p, smear=sm, z0=GRAPHITE_BE * g_log,
                                                   free=lambda z: gl(e, z, width(GRAPHITE_W, p)),
                                                   hint=hints.get(label), kq=kq))
            hints[label] = (got[1][0], got[1][3])
            return got
        return pick(label, lambda p: fit_region("C1s", e, y, st, p, smear=sm, kq=kq,
                                                lin=[gl(e, GRAPHITE_BE, width(GRAPHITE_W, p))]))

    def smear_of(region, par):
        """Offsets of the build-up charge from its region mean at eight instants of the region."""
        if par is None:
            return (0.0,)
        t0, t1 = span[region]
        q = par[1] + (par[0] - par[1]) * np.exp(-(t0 + (np.arange(8) + 0.5) * (t1 - t0) / 8) / par[2])
        return tuple(q - q.mean()) if np.ptp(q) > 0.02 else (0.0,)

    def broad_of(region, par, kap):
        """Uneven-charging FWHM, 2.3548 kappa |q|, at the instants smear_of uses (the region mean when it has none)."""
        if par is None or kap <= 0:
            return None
        t0, t1 = span[region]
        q = par[1] + (par[0] - par[1]) * np.exp(-(t0 + (np.arange(8) + 0.5) * (t1 - t0) / 8) / par[2])
        if np.ptp(q) <= 0.02:
            q = np.array([q.mean()])
        return tuple(2.3548 * kap * np.abs(q))

    par, hints, kappa = None, {}, 0.0
    rounds = 3 if uneven else (2 if drift == "build-up" else 1)
    for it in range(rounds):
        if par is not None and par[0] == par[1] and kappa <= 0.0:
            break                       # steady charge and even charging: nothing left to refine
        kq = {lab: broad_of(lab, par, kappa) for lab in span}
        if substrate:                   # graphite position free on the raw axis: pass energy and scale
            zs = [c1s(lab, None, smear_of(lab, par), kq[lab])[1][3] for lab in ("C1s", "C1s_end")]
        # the grounded graphite sits at (1 + g) 284.40 eV on the axis
        scale = float(np.mean(zs)) / GRAPHITE_BE if gain == "graphite" and substrate else g_log
        (pe_c, (pa, c_area, _, _)), (pe_cb, (pb, c_area_end, _, _)) = (c1s(lab, scale, smear_of(lab, par), kq[lab])
                                                          for lab in ("C1s", "C1s_end"))
        qa, qb = pa - 284.8, pb - 284.8
        eo, yo = arr["O1s"][:, 0] / scale, arr["O1s"][:, 1]
        sm_o = smear_of("O1s", par)
        pe_o, (o_pos, o_area, _, _) = pick("O1s", lambda p: fit_region("O1s", eo, yo, STEP / scale, p, smear=sm_o,
                                                                       kq=kq["O1s"]))
        qo = o_pos - O_LATTICE
        if drift == "build-up":
            if max(abs(qa - qb), abs(qa - qo)) < 0.03:
                par = (np.mean([qa, qb, qo]),) * 2 + (600.0,)
            else:
                obs = [("C1s", qa), ("O1s", qo), ("C1s_end", qb)]
                res = lambda p: [mean_charge(p[0], p[1], p[2], *span[k]) - q for k, q in obs]
                par = tuple(least_squares(res, [qa, qb, 400.0], bounds=([-5, -5, 20.0], [10, 10, 2e4])).x)
        else:
            par = (qa, qa, 600.0)
        if uneven and it < rounds - 1:
            # kappa from the two C 1s regions, whose pass energy the uncharged graphite fixes
            def cost(kap):
                tot = 0.0
                for lab in ("C1s", "C1s_end"):
                    e_, y_ = arr[lab][:, 0] / scale, arr[lab][:, 1]
                    p_ = chosen[lab]
                    lin = [gl(e_, GRAPHITE_BE, width(GRAPHITE_W, p_))] if substrate else ()
                    tot += fit_region("C1s", e_, y_, STEP / scale, p_, smear=smear_of(lab, par),
                                      kq=broad_of(lab, par, kap), lin=lin, hint=pa if lab == "C1s" else pb)[2]
                return tot
            kappa = float(minimize_scalar(cost, bounds=(0.0, 0.4), method="bounded", options={"xatol": 2e-3}).x)
            chosen.pop("O1s", None)     # chosen again with the broadening
    if drift == "build-up":
        def charge_at(region):
            return mean_charge(*par, *span[region])
    elif drift == "linear":
        mid = {k: 0.5 * sum(v) for k, v in span.items()}

        def charge_at(region):
            return qa + (qb - qa) * (mid[region] - mid["C1s"]) / (mid["C1s_end"] - mid["C1s"])
    else:
        def charge_at(region):
            return qa

    def tpe(p):
        return p if own_pass else 50

    e, y = arr["Ni2p3"][:, 0] / scale, arr["Ni2p3"][:, 1]
    q_ni = charge_at("Ni2p3")
    sm_ni = smear_of("Ni2p3", par if drift == "build-up" else None)
    kq_ni = broad_of("Ni2p3", par, kappa)
    pe_ni, (pos, area, _, _) = pick("Ni2p3", lambda p: fit_region("Ni2p3", e, y, STEP / scale, p, smear=sm_ni,
                                                                   kq=kq_ni))
    ni_be = pos - q_ni
    ni_free = area / dur["Ni2p3"] / transmission(tpe(pe_ni), HV - ni_be)

    e, y = arr["Fe2p3"][:, 0] / scale, arr["Fe2p3"][:, 1]
    q_fe = charge_at("Fe2p3")
    sm_fe = smear_of("Fe2p3", par if drift == "build-up" else None)
    kq_fe = broad_of("Fe2p3", par, kappa)
    kq_aug = kq_fe if kq_fe is not None and len(kq_fe) == len(sm_fe) else (
        (float(np.mean(kq_fe)),) * len(sm_fe) if kq_fe is not None else (0.0,) * len(sm_fe))

    def fe_fit(p):
        extra = None
        if auger != "none":
            pe_aug = tpe(p) if auger == "own" else 50
            aug_area = AUGER_RATIO * ni_free * transmission(pe_aug, HV - AUGER_BE) * dur["Fe2p3"]
            prof = sum(gl(e, AUGER_BE + q_fe + dq, float(np.hypot(width(AUGER_W, p), b)))
                       for dq, b in zip(sm_fe, kq_aug)) / len(sm_fe)
            extra = prof * aug_area / (prof.sum() * STEP / scale)
        return fit_region("Fe2p3", e, y, STEP / scale, p, charge=q_fe, extra=extra, smear=sm_fe, kq=kq_aug if kq_fe else None)
    pe_fe, (fe_pos, fe_area, _, _) = pick("Fe2p3", fe_fit)
    fe_ke = HV - (fe_pos - q_fe)
    n_ni = ni_free / RSF["Ni2p3"]
    n_fe = fe_area / dur["Fe2p3"] / transmission(tpe(pe_fe), fe_ke) / RSF["Fe2p3"]
    if overlayer:
        n_o = o_area / dur["O1s"] / transmission(tpe(pe_o), HV - O_LATTICE) / RSF["O1s"]
        n_c = c_area / dur["C1s"] / transmission(tpe(pe_c), HV - 284.8) / RSF["C1s"]
        ke = [HV - ni_be, fe_ke, HV - O_LATTICE]
        d = overlayer_thickness(n_c, [n_ni, n_fe, n_o], HV - 284.8, ke)
        if growth:
            # the overlayer may grow under the beam: the opening and closing C 1s give its thickness at both ends,
            # linear in time between them, and every oxide region is corrected at its own mid-time
            n_cb = c_area_end / dur["C1s_end"] / transmission(tpe(pe_cb), HV - 284.8) / RSF["C1s"]
            mids = {k: 0.5 * (v[0] + v[1]) for k, v in span.items()}
            ta, tb = mids["C1s"], mids["C1s_end"]
            t_ox = [mids["Ni2p3"], mids["Fe2p3"], mids["O1s"]]
            lam_c = imfp(HV - 284.8)

            def ends(dd):
                at = lambda tt: dd[0] + (dd[1] - dd[0]) * (tt - ta) / (tb - ta)
                tot = sum(n * np.exp(at(tt) / imfp(k)) for n, k, tt in zip([n_ni, n_fe, n_o], ke, t_ox))
                return [RHO_C * (1 - np.exp(-dd[0] / lam_c)) * tot - n_c, RHO_C * (1 - np.exp(-dd[1] / lam_c)) * tot - n_cb]
            da, db = least_squares(ends, [d, d], bounds=([0.0, 0.0], [20.0, 20.0])).x
            d_at = lambda tt: da + (db - da) * (tt - ta) / (tb - ta)
            n_ni, n_fe = n_ni * np.exp(d_at(t_ox[0]) / imfp(ke[0])), n_fe * np.exp(d_at(t_ox[1]) / imfp(ke[1]))
        else:
            n_ni, n_fe = n_ni * np.exp(d / imfp(ke[0])), n_fe * np.exp(d / imfp(ke[1]))
    return {"ni2p3_binding_ev": float(ni_be), "fe_fraction": float(n_fe / (n_fe + n_ni))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_electrode(sys.argv[1])))
