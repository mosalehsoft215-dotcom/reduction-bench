"""Reduces one lead isotope-dilution determination to the sample's lead mass fraction (dry mass) and its 206Pb/204Pb.

    python3 idms_reduce.py <determination directory> <result.json>

1. Count rates.  Every reading's rates are taken back to the counter's observed rates with the dead time the software
   applied (N_obs = N_rec / (1 + N_rec tau_set), tau_set from the sheet) and corrected with the detector's dead time
   tau (N = N_obs / (1 - N_obs tau)).
2. Zero of a solution.  The session's first rinse holds the baseline alone.  Memory washes out exponentially with one
   time constant T for all isotopes, so the two rinses read 0.75 min apart before a solution give T (the median over
   the session of 0.75 / ln(m1 / m2), m1 and m2 the two rinses less the baseline, on 206 and 208 where m2 exceeds
   300 cps) and the solution's zero is the baseline plus m2 carried forward to the solution's reading window (the mean
   of the exponential over that window against its mean over the second rinse's).
3. Dead time and mass bias together.  The assay standard is read at a low and a high level in turn; with the right
   dead time the 207/206 and 208/206 of all its readings follow one exponential-law mass bias f(t) = a + b t + c t^2
   (t the minute).  For each trial tau the linear least-squares fit of ln(R_meas / R_cert) = f(t) ln(m_i / m_206)
   (weights the counts of the smaller isotope) leaves a sum of squares; tau (0 to 100 ns) minimises it, and f(t) is
   the fit at that tau.
4. Mercury.  204Hg = net 202 x 0.2293 x (m_204 / m_202)^f taken off mass 204; ratios corrected by (m_i / m_206)^f.
5. Spike.  The working spike may hold laboratory lead (the assay standard's composition) besides the certified spike
   lead: atom fractions y = (1 - s) y_cert + s z, amount content c_Y.  A reverse blend of m_Z of the assay standard
   (c_Z) with m_Y of spike is (m_Z c_Z + s m_Y c_Y) of lead z plus (1 - s) m_Y c_Y of lead y_cert, so each of its
   ratios R (207/206, 208/206) gives rho = (y_cert,i - R y_cert,6) / (R z_6 - z_i) = (m_Z c_Z / m_Y) u + v with
   u = 1 / ((1 - s) c_Y) and v = s / (1 - s); the two reverse blends (a spike-poor and a spike-rich one) give u and v
   by least squares over both ratios, hence s and c_Y.
6. Blank: the spiked blank (laboratory lead) N_blank = N_Y (y_8 - R y_6) / (R z_6 - z_8), N_Y its diluted spike.
7. Sample: each blend's N_x = [N_Y (y_8 - R y_6) - N_blank (R z_6 - z_8)] / (R a_6 - a_8); the unspiked digest's
   fractions less its own blank give a (iterated with the mass fraction); mass fraction N_x M_x / (m_x (1 - water)),
   the water from the moisture portion, averaged over the three blends; 206Pb/204Pb = a_206 / a_204.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

MASS = {"202": 201.970643, "204": 203.973044, "206": 205.974465, "207": 206.975897, "208": 207.976652}
ISO = ("204", "206", "207", "208")
HG = 0.2293


def fractions(r):
    tot = sum(float(r[i]) for i in ISO)
    return {i: float(r[i]) / tot for i in ISO}


def lnm(i):
    return np.log(MASS[i] / MASS["206"])


def reduce_determination(folder, deadtime="fit", bias="drift", zero="washout", spike="solve", mercury=True,
                         moisture=True, blank=True, vessel=True, vessel_nat=True):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py: deadtime='sheet'
    keeps the software's dead time, bias='morning' takes the sheet's f for the whole session, zero='last' takes the
    last rinse before a solution as its zero, spike='certificate' takes the spike's certified ratios (its amount
    content the mean of the two reverse blends), spike='label' also its labelled amount content, mercury=False leaves
    204Hg in, moisture=False reports on the mass as weighed, blank=False leaves the blank in, vessel=False takes the
    mean of the three sample blends whatever their agreement, vessel_nat=False takes the unspiked digest's lead as the
    sample's and the procedural blank's only."""
    folder = Path(folder)
    sheet = json.loads((folder / "bench_sheet.json").read_text())
    with open(folder / "session.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    tau_set = sheet["instrument"]["dead_time_ns"] * 1e-9
    keys = ("202",) + ISO
    t_of = lambda r: float(r["minute"])  # noqa: E731
    dur = lambda r: float(r["seconds"]) / 60.0  # noqa: E731

    def rates(row, tau):
        out = {}
        for k in keys:
            rec = float(row["cps_" + k])
            obs = rec / (1 + rec * tau_set)
            out[k] = obs / (1 - obs * tau)
        return out

    sols = [j for j, r in enumerate(rows) if r["solution"] != "rinse"]

    def zeros(tau):
        base = rates(rows[0], tau)
        mem = {}
        times = []
        for j in sols:
            r1, r2 = rates(rows[j - 2], tau), rates(rows[j - 1], tau)
            m1 = {k: r1[k] - base[k] for k in keys}
            m2 = {k: r2[k] - base[k] for k in keys}
            mem[j] = m2
            for k in ("206", "208"):
                if m2[k] > 300 and m1[k] > m2[k]:
                    times.append((t_of(rows[j - 1]) - t_of(rows[j - 2])) / np.log(m1[k] / m2[k]))
        T = float(np.median(times)) if times else 1.0
        out = {}
        for j in sols:
            r2, sol = rows[j - 1], rows[j]
            if zero == "last":
                out[j] = {k: base[k] + mem[j][k] for k in keys}
                continue
            a2, b2 = t_of(r2), t_of(r2) + dur(r2)
            a, b = t_of(sol), t_of(sol) + dur(sol)
            carry = ((np.exp(-(a - a2) / T) - np.exp(-(b - a2) / T)) / (b - a)) / (
                (1 - np.exp(-(b2 - a2) / T)) / (b2 - a2))
            out[j] = {k: base[k] + mem[j][k] * carry for k in keys}
        return out

    def nets(tau):
        z0 = zeros(tau)
        out = {}
        for j in sols:
            r = rates(rows[j], tau)
            out[j] = {k: r[k] - z0[j][k] for k in keys}
        return out

    cert = sheet["assay_standard"]["ratios_to_206"]
    assays = [j for j in sols if rows[j]["solution"] == "assay"]

    def bias_fit(tau):
        n = nets(tau)
        A, yv, w = [], [], []
        for j in assays:
            tt = t_of(rows[j]) + dur(rows[j]) / 2
            for i in ("207", "208"):
                A.append(np.array([1.0, tt, tt * tt]) * lnm(i))
                yv.append(np.log(n[j][i] / n[j]["206"] / cert[i]))
                w.append(np.sqrt(min(n[j][i], n[j]["206"]) * float(rows[j]["seconds"])))
        A, yv, w = np.array(A), np.array(yv), np.array(w)
        coef, *_ = np.linalg.lstsq(A * w[:, None], yv * w, rcond=None)
        return float(np.sum(((A @ coef - yv) * w) ** 2)), coef

    if deadtime == "sheet":
        tau = tau_set
    else:
        tau = minimize_scalar(lambda x: bias_fit(x * 1e-9)[0], bounds=(0.0, 100.0), method="bounded",
                              options={"xatol": 1e-4}).x * 1e-9
    coef = bias_fit(tau)[1]
    net = nets(tau)

    def f_at(j):
        if bias == "morning":
            return sheet["instrument"]["mass_bias_f_morning"]
        tt = t_of(rows[j]) + dur(rows[j]) / 2
        return float(coef[0] + coef[1] * tt + coef[2] * tt * tt)

    def ratios(label):
        j = next(k for k in sols if rows[k]["solution"] == label)
        n, f = net[j], f_at(j)
        n204 = n["204"] - (n["202"] * HG * (MASS["204"] / MASS["202"]) ** f if mercury else 0.0)
        out = {"204": n204 / n["206"], "206": 1.0, "207": n["207"] / n["206"], "208": n["208"] / n["206"]}
        return {i: out[i] / np.exp(f * lnm(i)) for i in ISO}

    z = fractions(cert)
    y_cert = fractions(sheet["spike"]["ratios_to_206"])
    weigh = {w["blend"]: w for w in sheet["weighings"]}
    c_z = sheet["assay_standard"]["nmol_per_g"]
    revs = sorted(k for k in weigh if k.startswith("reverse-"))
    if spike == "solve":
        lhs, rhs = [], []
        for lab in revs:
            r = ratios(lab)
            for i in ("207", "208"):
                rhs.append((y_cert[i] - r[i] * y_cert["206"]) / (r[i] * z["206"] - z[i]))   # lead z per lead y_cert
                lhs.append([weigh[lab]["assay_g"] * c_z / weigh[lab]["spike_g"], 1.0])
        u, v = np.linalg.lstsq(np.array(lhs), np.array(rhs), rcond=None)[0]
        s = v / (1 + v)
        c_y = 1.0 / ((1 - s) * u)
        y = {i: (1 - s) * y_cert[i] + s * z[i] for i in ISO}
    else:
        y = y_cert
        if spike == "label":
            c_y = sheet["spike"]["label_nmol_per_g"]
        else:
            cs = []
            for lab in revs:
                rr = ratios(lab)["208"]
                n_y = weigh[lab]["assay_g"] * c_z * (z["208"] - rr * z["206"]) / (rr * y["206"] - y["208"])
                cs.append(n_y / weigh[lab]["spike_g"])
            c_y = float(np.mean(cs))
    n_blank = 0.0
    if blank:
        rb = ratios("blank")["208"]
        n_yb = weigh["blank"]["spike_g"] * c_y / weigh["blank"]["spike_dilution"]
        n_blank = n_yb * (y["208"] - rb * y["206"]) / (rb * z["206"] - z["208"])
    mo = sheet["moisture"]
    water = 1 - mo["after_g"] / mo["before_g"] if moisture else 0.0
    mix = fractions(ratios("natural"))
    blends = {lab: ratios(lab)["208"] for lab in weigh if lab.startswith("blend-")}
    r207 = {lab: ratios(lab)["207"] for lab in blends}

    def amount(lab, a, i, rr):
        n_y = weigh[lab]["spike_g"] * c_y
        return (n_y * (y[i] - rr * y["206"]) - n_blank * (rr * z["206"] - z[i])) / (rr * a["206"] - a[i])

    def composition(extra, n_nat):          # the sample's fractions from the unspiked digest, extra vessel lead out
        return {i: ((n_nat + n_blank + extra) * mix[i] - (n_blank + extra) * z[i]) / n_nat for i in ISO}

    a, w_dry, kept, extra = dict(mix), None, sorted(blends), 0.0
    for _ in range(6):
        m_x = sum(a[i] * MASS[i] for i in ISO)
        out = []
        for lab, rr in sorted(blends.items()):
            n_y = weigh[lab]["spike_g"] * c_y
            n_x = (n_y * (y["208"] - rr * y["206"]) - n_blank * (rr * z["206"] - z["208"])) / (rr * a["206"] - a["208"])
            out.append(n_x * m_x * 1e-3 / (weigh[lab]["sample_g"] * (1 - water)))
        out = np.array(out)
        kept = sorted(blends)
        if vessel and out.size == 3:                 # a digest whose vessel added lead stands apart from the others
            dev = np.abs(out / np.median(out) - 1)
            far = int(np.argmax(dev))
            pair = np.delete(out, far)
            if dev[far] > 0.01 and abs(pair[0] / pair[1] - 1) < 0.003:
                out = pair
                kept = [lab for k, lab in enumerate(sorted(blends)) if k != far]
        w_dry = float(np.mean(out))
        n_nat = weigh["natural"]["sample_g"] * (1 - water) * w_dry * 1e3 / m_x
        if vessel_nat:   # lead the unspiked digest's vessel added: the one that makes 207 and 208 tell the same amount
            def disagree(e):
                aa = composition(e, n_nat)
                return float(sum((amount(lab, aa, "207", r207[lab]) / amount(lab, aa, "208", blends[lab]) - 1) ** 2
                                 for lab in kept))
            extra = minimize_scalar(disagree, bounds=(0.0, 0.4 * n_nat), method="bounded",
                                    options={"xatol": 1e-6 * n_nat}).x
        a = composition(extra, n_nat)
    return {"lead_ug_per_g": w_dry, "pb206_pb204": float(a["206"] / a["204"])}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_determination(sys.argv[1])))
