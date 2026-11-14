"""Reference fit of the separator field trial: permeances of CO2 and CH4 at 35 C (and of C2H6 and N2, with all four
activation energies), from the raw plant records.

    python permeance_fit.py <data folder> <results.json>

Records: historian.csv (DCS, 10 s), gc_analyses.csv (one stream per analysis), operator_log.txt.
GC: the CO2 response is corrected with the span-gas analyses (stream 4, straight line in time through all of them),
then every analysis is normalised again; each analysis is dated back to the moment its gas left the process
(end of analysis less the analysis time less the transport time of its sample line).
FT-103: zero offset taken as a straight line between the two zero checks of the operator log; the reading is then
multiplied by the permeate's gas correction factor 1 / sum(y_i / K_i). FT-101 and FT-102 are standard m3/h at 15 C.
Steady data: each operating point after a settling time, split where the feed supply changes wells (a switch and its
blend are left out).
Fibre break: when the operator heard a whistling noise (or, if nobody did and the inspection found broken fibres, at
the largest step of FT-103 inside a settled operating point), FT-103 jumps at constant conditions; that fixes the
moment. From then on the broken fibres no longer separate, and shell gas of the composition present where they broke
leaves the shell there through their open bores into the permeate header: laminar, compressible, at the local gas
temperature along the bore and in the potted bores of the open end,
n = N_b pi d_i^4 (p_feed^2 - p_permeate^2) / (256 mu R integral(T dz)), mu by Wilke's rule for that gas.
For every steady piece the historian is averaged around the permeate sample times; feed composition is the mean of the
feed analyses of the same well period, retentate and permeate compositions the means of their analyses in the piece.
Separator: feed outside the fibres from the feed end; the permeate inside the bores leaves at the end where the bores
are open (campaign.json flow_paths: the feed end for counter-current modules, the retentate end for co-current ones).
Area: outer surface of the open fibres (built less plugged) over the length between the pottings. Bore pressure: rises
from the nozzle through the potted bores of the open end (gas temperature of that end) and along the fibres towards
the closed ends (laminar, Wilke viscosity of the mixed permeate present there). Temperature: linear between TT-1 and the separator outlet,
which is TT-2 corrected for the throttling across PCV-102 (same enthalpy, Peng-Robinson departure + ideal-gas heat
capacity). Driving force: Peng-Robinson fugacity difference.
Fit: weighted least squares of the retentate and permeate component flows of all steady pieces.
"""
import csv
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

R = 8.314462618
GPU = 3.348e-10
T_REF = 308.15
TC = np.array([304.13, 190.56, 305.32, 126.19])
PC = np.array([73.773e5, 45.992e5, 48.722e5, 33.958e5])
OMEGA = np.array([0.2239, 0.0115, 0.0995, 0.0372])
KIJ = np.array([[0.0, 0.0919, 0.1322, -0.0170], [0.0919, 0.0, 0.0026, 0.0311],
                [0.1322, 0.0026, 0.0, 0.0515], [-0.0170, 0.0311, 0.0515, 0.0]])
# module_datasheet.txt
FIBRES, OD, ID = 33870 - 2070, 300e-6, 150e-6       # set per trial by configure() from campaign.json
PLUGGED = 2070
TOTAL_LENGTH, POTTING = 1.340, 0.120
# test_instruments.txt
K_THERMAL = np.array([0.74, 0.72, 0.50, 1.00])
V_STD, V_NORM = R * 288.15 / 101325.0, R * 273.15 / 101325.0
# viscosities for Wilke's rule (Pa s at 35 C) and molar masses
MU = np.array([1.53e-5, 1.13e-5, 0.95e-5, 1.80e-5])      # peng_robinson.csv, viscosity column
MW = np.array([44.01, 16.04, 30.07, 28.01])

DEFAULT = dict(area="outer", length="active", bore_drop=True, potting_drop=True, valve_jt=True, permeate="gcf", driving="fugacity", n=60, sealed=True,
               temperature="profile", gauge=True, gc_span=True, zero=True, transport=True, wells=True, settle=60.0, fibre_break="correct",
               force=None, find_break=True, leak="local")
# gc_analyses.csv and test_instruments.txt
GC_CYCLE = 4.0                                  # min from injection to the end of an analysis
TRANSPORT = {1: 1.0, 2: 2.4, 3: 6.8}            # min, sample lines: feed, retentate, permeate
SPAN_GAS = np.array([10.00, 84.00, 4.00, 2.00])  # mol %, certificate
BLEND, MARGIN = 15.0, 3.0                        # min


def fugacity(y, p, t, driving="fugacity"):
    p = np.asarray(p, float)[..., None]
    if driving != "fugacity":
        return y * p
    t = np.asarray(t, float)[..., None]
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA ** 2
    a = 0.45724 * R ** 2 * TC ** 2 / PC * (1 + kappa * (1 - np.sqrt(t / TC))) ** 2
    b = 0.07780 * R * TC / PC
    aij = np.sqrt(a[..., :, None] * a[..., None, :]) * (1 - KIJ)
    sa = np.einsum("...ij,...j->...i", aij, y)
    am, bm = (y * sa).sum(-1, keepdims=True), (y * b).sum(-1, keepdims=True)
    A, B = am * p / (R * t) ** 2, bm * p / (R * t)
    z = np.ones_like(A)
    for _ in range(12):
        f = z ** 3 - (1 - B) * z ** 2 + (A - 3 * B ** 2 - 2 * B) * z - (A * B - B ** 2 - B ** 3)
        z = z - f / (3 * z ** 2 - 2 * (1 - B) * z + (A - 3 * B ** 2 - 2 * B))
    s2 = np.sqrt(2.0)
    lnphi = b / bm * (z - 1) - np.log(z - B) - A / (2 * s2 * B) * (2 * sa / am - b / bm) * np.log((z + (1 + s2) * B) / (z + (1 - s2) * B))
    return y * p * np.exp(lnphi)


# peng_robinson.csv: ideal-gas heat capacity Cp = a + bT + cT^2 + dT^3 (J/(mol K))
CP = np.array([[19.80, 7.344e-2, -5.602e-5, 1.715e-8], [19.25, 5.213e-2, 1.197e-5, -1.132e-8],
               [5.409, 1.781e-1, -6.938e-5, 8.713e-9], [31.15, -1.357e-2, 2.680e-5, -1.168e-8]])


def molar_enthalpy(y, p, t):
    """ideal-gas enthalpy from 298.15 K plus the Peng-Robinson departure, J/mol"""
    t0 = 298.15
    powers = np.array([t - t0, (t ** 2 - t0 ** 2) / 2, (t ** 3 - t0 ** 3) / 3, (t ** 4 - t0 ** 4) / 4])
    ideal = y @ (CP @ powers)
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA ** 2
    s = 1 + kappa * (1 - np.sqrt(t / TC))
    ac = 0.45724 * R ** 2 * TC ** 2 / PC
    a, da = ac * s ** 2, -ac * kappa * s / np.sqrt(t * TC)
    g = np.sqrt(np.outer(a, a))
    am = y @ (g * (1 - KIJ)) @ y
    dam = y @ ((np.outer(da, a) + np.outer(a, da)) / (2 * g) * (1 - KIJ)) @ y
    bm = y @ (0.07780 * R * TC / PC)
    A, Bv = am * p / (R * t) ** 2, bm * p / (R * t)
    roots = np.roots([1, -(1 - Bv), A - 3 * Bv ** 2 - 2 * Bv, -(A * Bv - Bv ** 2 - Bv ** 3)])
    z = max(r.real for r in roots if abs(r.imag) < 1e-9)
    s2 = np.sqrt(2.0)
    return ideal + R * t * (z - 1) + (t * dam - am) / (2 * s2 * bm) * np.log((z + (1 + s2) * Bv) / (z + (1 - s2) * Bv))


def before_valve(y, p_up, p_down, t_down):
    """temperature upstream of a throttling valve, from the temperature measured downstream (same enthalpy)"""
    h = molar_enthalpy(y, p_down, t_down)
    lo, hi = t_down - 5.0, t_down + 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if molar_enthalpy(y, p_up, mid) < h else (lo, mid)
    return 0.5 * (lo + hi)


def wilke(y):
    """viscosity of a gas mixture (Pa s), y (..., 4)"""
    phi = (1 + np.sqrt(MU[:, None] / MU[None, :]) * (MW[None, :] / MW[:, None]) ** 0.25) ** 2 / np.sqrt(8 * (1 + MW[:, None] / MW[None, :]))
    return (y * MU / np.einsum("...j,ij->...i", y, phi)).sum(-1)


def module(runs, q35, e, o):
    """retentate (shell outflow) and permeate (bores plus leak) component flows (mol/s) of all runs"""
    feed = runs["feed"]
    nr = len(feed)
    n = o["n"]
    length = TOTAL_LENGTH - 2 * POTTING if o["length"] == "active" else TOTAL_LENGTH
    outlet = runs["outlet"]
    count = runs["count"]                          # broken fibres in each run (0 before the break)
    zf = np.linspace(0.0, 1.0, n + 1)
    kb = None
    if count.max() > 0:
        sb = 0.0 if o["leak"] == "feed_gas" else runs["zb"] / length
        zf = np.unique(np.round(np.r_[zf, sb], 12))
        kb = int(np.argmin(np.abs(zf - sb)))
    nz = zf.size
    dz = np.diff(zf) * length
    fibres = (FIBRES if o["sealed"] else FIBRES + PLUGGED) - count
    da = (fibres * np.pi * (OD if o["area"] == "outer" else ID))[:, None]
    hp = 256 * R / (np.pi * ID ** 4 * fibres)
    if o["temperature"] == "profile":       # linear between the inlet and outlet thermocouples
        t = runs["t"][:, None] + (runs["t_out"] - runs["t"])[:, None] * zf
    elif o["temperature"] == "mean":
        t = (0.5 * (runs["t"] + runs["t_out"]))[:, None] * np.ones(nz)
    else:
        t = runs["t"][:, None] * np.ones(nz)
    q = q35 * GPU * np.exp(-e * 1e3 / R * (1 / t[..., None] - 1 / T_REF))
    pf = runs["p_feed"][:, None] * np.ones(nz)
    yb = np.repeat((feed / feed.sum(1, keepdims=True))[:, None, :], nz, 1)
    pb = runs["p_perm"][:, None] * np.ones(nz)
    g = np.zeros_like(yb)
    lk = np.zeros_like(feed)
    if kb is not None:      # integral of T along the open bore from the break to the nozzle
        tf, tr = t[:, 0], t[:, -1]
        if o["leak"] == "feed_gas":         # misreading: feed gas at the inlet temperature, taken off before the membrane
            tint = (runs["zb"] + POTTING) * tf
        elif outlet == "feed":
            tint = length * (sb * tf + (tr - tf) * sb ** 2 / 2) + POTTING * tf
        else:
            tint = length * ((1 - sb) * tf + (tr - tf) * (1 - sb ** 2) / 2) + POTTING * tr
        cl = count * np.pi * ID ** 4 * (runs["p_feed"] ** 2 - runs["p_perm"] ** 2) / (256 * R * tint)
    for sweep in range(300):
        fb = fugacity(yb, pb, t, o["driving"])
        flow = np.zeros_like(yb)
        flow[:, 0] = feed
        flux = np.zeros_like(yb)
        for k in range(nz):
            x = flow[:, k] / flow[:, k].sum(1, keepdims=True)
            flux[:, k] = q[:, k] * (fugacity(x, pf[:, k], t[:, k], o["driving"]) - fb[:, k]) * da
            if k == kb:     # shell gas of the local composition leaves through the broken bores
                lk = (cl / wilke(x))[:, None] * x
                flow[:, k] = flow[:, k] - lk
            if k < nz - 1:
                pred = flow[:, k] - flux[:, k] * dz[k]
                xp = pred / pred.sum(1, keepdims=True)
                nxt = q[:, k + 1] * (fugacity(xp, pf[:, k + 1], t[:, k + 1], o["driving"]) - fb[:, k + 1]) * da
                flow[:, k + 1] = flow[:, k] - 0.5 * (flux[:, k] + nxt) * dz[k]
        seg = 0.5 * (flux[:, 1:] + flux[:, :-1]) * dz[None, :, None]
        g_new = np.zeros_like(g)
        if outlet == "feed":        # bores closed at the retentate end
            g_new[:, :-1] = np.cumsum(seg[:, ::-1], 1)[:, ::-1]
        else:                       # bores closed at the feed end
            g_new[:, 1:] = np.cumsum(seg, 1)
        tot = g_new.sum(-1, keepdims=True)
        y_new = np.where(tot > 0, g_new / np.where(tot > 0, tot, 1.0), flux / flux.sum(-1, keepdims=True))
        p_new = runs["p_perm"][:, None] * np.ones(nz)
        if o["bore_drop"]:
            mu = wilke(y_new)
            gm = 0.5 * (tot[:, 1:, 0] + tot[:, :-1, 0])
            mm = 0.5 * (mu[:, 1:] + mu[:, :-1])
            inc = hp[:, None] * mm * gm * 0.5 * (t[:, 1:] + t[:, :-1]) * dz[None, :]
            e_ = 0 if outlet == "feed" else -1
            p_edge2 = runs["p_perm"] ** 2
            if o["potting_drop"]:      # all of the permeate flows through the potted bores of the open end
                p_edge2 = p_edge2 + hp * mu[:, e_] * tot[:, e_, 0] * t[:, e_] * POTTING
            if outlet == "feed":
                p_new = np.sqrt(p_edge2[:, None] + np.concatenate([np.zeros((nr, 1)), np.cumsum(inc, 1)], 1))
            else:
                p_new = np.sqrt(p_edge2[:, None] + np.concatenate([np.cumsum(inc[:, ::-1], 1)[:, ::-1], np.zeros((nr, 1))], 1))
        change = np.abs(g_new - g).max() / g_new.max()
        g, yb, pb = g_new, 0.5 * yb + 0.5 * y_new, 0.5 * pb + 0.5 * p_new
        if change < 1e-11:
            break
    return flow[:, -1], (g[:, 0] if outlet == "feed" else g[:, -1]) + lk


def when(text):
    text = text.strip()
    return datetime.strptime(text, "%Y-%m-%d %H:%M:%S" if len(text) > 16 else "%Y-%m-%d %H:%M")


def records(folder):
    folder = Path(folder)
    with open(folder / "historian.csv") as fh:
        rows = list(csv.reader(fh))
    head, rows = rows[0], rows[1:]
    t0 = when(rows[0][0])
    minute = lambda d: (d - t0).total_seconds() / 60.0
    hist = {h: np.array([float(r[i]) if r[i] else np.nan for r in rows]) for i, h in enumerate(head) if i}
    hist["t"] = np.array([minute(when(r[0])) for r in rows])
    with open(folder / "gc_analyses.csv") as fh:
        gc = list(csv.DictReader(fh))
    gc = dict(end=np.array([minute(when(r["analysis_end"])) for r in gc]), stream=np.array([int(r["stream"]) for r in gc]),
              y=np.array([[float(r[f"{c}_molpct"]) for c in ("CO2", "CH4", "C2H6", "N2")] for r in gc]))
    log = []
    for line in (folder / "operator_log.txt").read_text().splitlines():
        if line.strip():
            log.append((minute(when(line[:16])), line[16:].strip()))
    return hist, gc, log


def step(t, x, near):
    """moment and size of the step in x near the time `near` (min): one quadratic trend plus a step, best split"""
    w = (t > near - 7.0) & (t < near + 5.0) & np.isfinite(x)
    tt, xx = t[w], x[w]
    best = None
    for j in range(6, len(tt) - 6):
        tb = 0.5 * (tt[j - 1] + tt[j])
        a = np.column_stack([np.ones_like(tt), tt - tb, (tt - tb) ** 2, (tt > tb).astype(float)])
        coef, sse = np.linalg.lstsq(a, xx, rcond=None)[:2]
        if best is None or sse[0] < best[0]:
            best = (sse[0], tb, coef[3])
    return best[1], best[2]


def find_step(t, x, sets, end):
    """time of the largest step of x inside the settled part of any operating point: 3 min means on either side"""
    best, when = 0.0, None
    for lo, hi in zip(sets, sets[1:] + [end]):
        idx = np.flatnonzero((t >= lo + 30.0) & (t <= hi - 5.0))
        for i in idx[18:-18]:
            if not (np.isfinite(x[i - 18:i]).any() and np.isfinite(x[i:i + 18]).any()):    # historian outage
                continue
            a, b = np.nanmean(x[i - 18:i]), np.nanmean(x[i:i + 18])
            if abs(b - a) / abs(a) > best:
                best, when = abs(b - a) / abs(a), t[i]
    return when


def configure(folder):
    """module geometry and instrument arrangement of this trial, from campaign.json"""
    global FIBRES, PLUGGED, OD, ID, TOTAL_LENGTH, POTTING
    cj = json.loads((Path(folder) / "campaign.json").read_text())
    m = cj["module"]
    PLUGGED = m["fibres_plugged_after_factory_leak_test"]
    FIBRES = m["fibres_built"] - PLUGGED
    OD, ID = m["fibre_od_um"] * 1e-6, m["fibre_id_um"] * 1e-6
    TOTAL_LENGTH, POTTING = m["fibre_length_between_outer_potting_faces_m"], m["potting_block_thickness_m"]
    ins, gcj = cj["instruments"], cj["gas_chromatograph"]
    gases = ("CO2", "CH4", "C2H6", "N2")
    found = cj["post_trial_inspection"]
    words = dict(one=1, two=2, three=3, four=4, five=5, six=6)
    m_break = re.search(r"(\w+) fibres? broken ([\d.]+) m from the (feed|retentate)-end", found)
    return dict(gauge_pt=ins["pressure_transmitters"].startswith("gauge"),
                tt2_after_valve="downstream of PCV-102" in ins["TT-2"],
                outlet="ret" if "closed in the feed-end potting" in m["flow_paths"] else "feed",
                break_end=({"feed": "feed", "retentate": "ret"}[m_break.group(3)] if m_break else "ret"),
                break_d=float(m_break.group(2)) if m_break else 0.02,
                break_count=(words.get(m_break.group(1), None) or int(m_break.group(1))) if m_break else 0,
                defect="broken" in found,
                span=[gcj["span_gas_certificate_mol_percent"][g] for g in gases],
                gcf=[ins["FT-103"]["gas_correction_factors"][g] for g in gases],
                transport={1: gcj["sample_line_transport_min"]["feed"], 2: gcj["sample_line_transport_min"]["retentate"],
                           3: gcj["sample_line_transport_min"]["permeate"]})


def read(folder, o):
    global SPAN_GAS, K_THERMAL, TRANSPORT
    arr = dict(configure(folder), **(o["force"] or {}))   # force: an arrangement taken against campaign.json
    SPAN_GAS, K_THERMAL, TRANSPORT = np.array(arr["span"], float), np.array(arr["gcf"], float), arr["transport"]
    hist, gc, log = records(folder)
    end = hist["t"][-1]
    sets = [m for m, s in log if s.startswith("Operating point")]
    switches = [m for m, s in log if s.startswith("Upstream supply changed")]
    zeros = [(m, float(s.split("reads ")[1].split()[0])) for m, s in log if "zero check" in s]
    broken = [m for m, s in log if "whistl" in s.lower()]
    if o["gauge"] and arr["gauge_pt"]:          # gauge transmitters: the DCS adds a standard atmosphere, the site barometer is lower
        baro = np.mean([float(s.split("barometer ")[1].split()[0]) / 1000.0 for m, s in log if "barometer" in s])
        for k in ("PT-101_bara", "PT-103_bara", "PT-104_bara"):
            hist[k] = hist[k] - 1.01325 + baro

    y = gc["y"].copy()
    if o["gc_span"]:        # CO2 response from the span gas, straight line in time
        sp = gc["stream"] == 4
        f = y[sp, 0] / y[sp, 1] * SPAN_GAS[1] / SPAN_GAS[0]
        y[:, 0] = y[:, 0] / np.polyval(np.polyfit(gc["end"][sp], f, 1), gc["end"])
    y = y / y.sum(1, keepdims=True)
    lag = np.array([GC_CYCLE + TRANSPORT.get(s, 0.0) for s in gc["stream"]]) if o["transport"] else 0.0
    taken = gc["end"] - lag

    zero = np.zeros_like(hist["t"])
    if o["zero"]:
        (m0, z0), (m1, z1) = zeros[0], zeros[-1]
        zero = z0 + (z1 - z0) * (hist["t"] - m0) / (m1 - m0)
    ft103 = hist["FT-103_Nm3h_N2"] - zero
    if arr["defect"] and not broken and o["find_break"]:    # fibres found broken, nobody noticed: the largest FT-103 step
        broken = [find_step(hist["t"], ft103, sets, end)]

    tb = None
    if o["fibre_break"] == "correct" and broken:     # the moment of the break: trend plus step of FT-103
        tb = step(hist["t"], ft103, broken[0])[0]

    last = len(sets)
    if o["fibre_break"] == "exclude" and broken:
        last = max(k for k, m in enumerate(sets) if m <= broken[0])
    bounds = sets + [end]
    edges = [-np.inf] + switches + [np.inf]
    wells = [(edges[i] + (BLEND + MARGIN if i else 0.0), edges[i + 1]) for i in range(len(edges) - 1)]
    pieces = []
    for k in range(last):
        lo, hi = bounds[k] + o["settle"], bounds[k + 1]
        if o["wells"]:
            pieces += [(max(lo, a), min(hi, b), a, b) for a, b in wells if min(hi, b) - max(lo, a) > 0]
        else:
            pieces.append((lo, hi, lo, hi))

    if tb is not None:      # no piece may straddle the break
        pieces = [q for lo, hi, wa, wb in pieces
                  for q in ([(lo, tb - 1.0, wa, wb), (tb + 2.0, hi, wa, wb)] if lo < tb < hi else [(lo, hi, wa, wb)])]
    keys = ("FT-101_Sm3h", "FT-102_Sm3h", "PT-101_bara", "PT-103_bara", "PT-104_bara", "TT-1_C", "TT-2_C")
    runs = []
    for lo, hi, wa, wb in pieces:
        pick = lambda s, a, b: (gc["stream"] == s) & (taken >= a) & (taken < b)
        ps, fs, rs = pick(3, lo, hi), pick(1, wa, wb), pick(2, lo, hi)
        if ps.sum() < 2 or not fs.any() or not rs.any():
            continue
        near = np.zeros_like(hist["t"], bool)
        for s in taken[ps]:
            near |= np.abs(hist["t"] - s) <= 1.0
        mean = {k: np.nanmean(hist[k][near]) for k in keys}
        mean["FT-103"] = np.nanmean(ft103[near])
        runs.append(dict(mean, feed_y=y[fs].mean(0), ret_y=y[rs].mean(0), perm_y=y[ps].mean(0), after=tb is not None and lo > tb))

    col = lambda k: np.array([r[k] for r in runs])
    fy, ry, py = col("feed_y"), col("ret_y"), col("perm_y")
    feed = col("FT-101_Sm3h") / 3600 / V_STD
    ret = col("FT-102_Sm3h") / 3600 / V_STD
    if o["permeate"] == "gcf":
        per = col("FT-103") / 3600 / V_NORM / (py / K_THERMAL).sum(1)
    elif o["permeate"] == "balance":
        per = feed - ret
    else:                                   # meter reading taken as the permeate flow
        per = col("FT-103") / 3600 / V_NORM
    t_tt2, p_feed, p_down = col("TT-2_C") + 273.15, col("PT-101_bara") * 1e5, col("PT-104_bara") * 1e5
    if o["valve_jt"] and arr["tt2_after_valve"]:       # TT-2 downstream of PCV-102: undo the isenthalpic expansion
        t_out = np.array([before_valve(ry[k], p_feed[k], p_down[k], t_tt2[k]) for k in range(len(runs))])
    else:
        t_out = t_tt2
    ret_m, per_m = ret[:, None] * ry, per[:, None] * py
    active = TOTAL_LENGTH - 2 * POTTING
    runs = dict(feed=feed[:, None] * fy, p_feed=p_feed, p_perm=col("PT-103_bara") * 1e5,
                t=col("TT-1_C") + 273.15, t_out=t_out, outlet=arr["outlet"],
                count=np.where(col("after"), arr["break_count"], 0) if tb is not None else np.zeros(len(p_feed), int),
                zb=arr["break_d"] if arr["break_end"] == "feed" else active - arr["break_d"])
    return runs, ret_m, per_m


def solve(folder, **change):
    o = dict(DEFAULT, **change)
    runs, ret_m, per_m = read(folder, o)
    scale = np.concatenate([np.maximum(ret_m, 1e-3 * ret_m.sum(1, keepdims=True)), np.maximum(per_m, 1e-3 * per_m.sum(1, keepdims=True))], 1)

    def residual(x):
        q35, e = np.exp(x[:4]), x[4:]
        ret, per = module(runs, q35, e, o)
        return ((np.concatenate([ret, per], 1) - np.concatenate([ret_m, per_m], 1)) / (0.01 * scale)).ravel()

    x0 = np.concatenate([np.log([50.0, 2.0, 0.6, 2.0]), [5.0, 15.0, 15.0, 10.0]])
    fit = least_squares(residual, x0, x_scale=np.r_[np.full(4, 0.1), np.full(4, 2.0)], xtol=1e-10, ftol=1e-10)
    q35, e = np.exp(fit.x[:4]), fit.x[4:]
    res = dict(Q_CO2_GPU=q35[0], Q_CH4_GPU=q35[1], E_CO2_kJ_mol=e[0], E_CH4_kJ_mol=e[1])
    info = dict(pieces=len(runs["p_feed"]), Q_GPU=q35.tolist(), E_kJ_mol=e.tolist(), rms_residual_percent=float(np.sqrt(np.mean(fit.fun ** 2))), evaluations=int(fit.nfev))
    return {k: float(v) for k, v in res.items()}, info


if __name__ == "__main__":
    results, details = solve(sys.argv[1])
    Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
    Path(sys.argv[2]).write_text(json.dumps(results, indent=1) + "\n")
    print(json.dumps(results, indent=1))
    print(json.dumps(details, indent=1))
