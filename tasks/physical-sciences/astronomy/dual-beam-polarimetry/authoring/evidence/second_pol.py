"""Independent second reduction, written apart from solution/pol_reduce.py, for the allowances.

    python3 second_pol.py <visit folder> <output file>

Where it differs from the reference:
- each beam's flux is taken to infinite aperture: a Moffat curve of growth is fitted to the aperture sums that no
  neighbouring image reaches (when at least three are clean; otherwise the largest clean one is used as it is);
- the position angle of each target-field exposure is the least-squares rotation that carries the stars'
  detector vectors onto their catalogue sky vectors;
- the Stokes fit is nonlinear, on log(o / e) = log t_j + log((1 + m) / (1 - m)) with scipy least squares;
- the efficiency is a free parameter of the polarised standard's fit, with the Serkowski degree in the band;
- the unpolarised standard's two rotator settings are fitted with an instrument-fixed and a sky-fixed vector;
- off axis it works in two stages: each field star's own apparent instrument-frame Stokes first, then a
  regression of those on a common vector plus a radial template rotated into the instrument frame.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def growth_total(radii, sums):
    """Total flux from aperture sums: Moffat curve of growth with free width and index, fitted in flux."""
    radii = np.asarray(radii, float)
    sums = np.asarray(sums, float)

    def model(p):
        tot, a, beta = p
        return tot * (1.0 - (1.0 + (radii / a) ** 2) ** (1.0 - beta))

    start = [sums[-1] * 1.01, radii[1], 4.0]
    fit = least_squares(lambda p: (model(p) - sums) / np.sqrt(np.abs(sums) + 1.0), start,
                        bounds=([0.0, 0.3, 1.5], [np.inf, 60.0, 12.0]))
    return fit.x[0]


def serkowski(p_max, lam_max, k, lam):
    return p_max * np.exp(-k * np.log(lam_max / lam) ** 2)


def records(folder):
    sheet = json.loads((folder / "visit_sheet.json").read_text())
    with open(folder / "exposures.csv") as f:
        expo = {r["exposure"]: r for r in csv.DictReader(f)}
    with open(folder / "photometry.csv") as f:
        phot = list(csv.DictReader(f))
    with open(folder / "field_catalogue.csv") as f:
        cat = {r["source"]: np.array([float(r["east_arcsec"]), float(r["north_arcsec"])]) for r in csv.DictReader(f)}
    radii = sorted(sheet["aperture_radii_px"])
    phi0 = sheet["hwp_zero_offset_deg"][sheet["filter"]]
    x0, y0 = sheet["optical_axis_px"]
    sx, sy = sheet["e_image_offset_px"]
    fw = max(float(e["seeing_arcsec"]) for e in expo.values()) / sheet["pixel_scale_arcsec"] * 1.15
    # clean radii per source: no other o or e image closer than r + 1.6 seeing widths to either own image
    pos = {}
    for r in phot:
        pos.setdefault((expo[r["exposure"]]["object"], r["source"]), []).append((float(r["x_px"]), float(r["y_px"])))
    pos = {k: np.median(np.array(v), axis=0) for k, v in pos.items()}
    clean = {}
    for (obj, src), c in pos.items():
        own = np.array([c, c + [sx, sy]])
        oth = np.array([q for (o2, s2), c2 in pos.items() if o2 == obj and s2 != src for q in (c2, c2 + [sx, sy])] or [[1e9, 1e9]])
        gap = np.min(np.hypot(own[:, None, 0] - oth[None, :, 0], own[:, None, 1] - oth[None, :, 1]))
        clean[(obj, src)] = [k for k in radii if k + 1.6 * fw <= gap] or radii[:1]
    # position angle of each target-field exposure: least-squares rotation of detector vectors onto sky vectors
    angle = {}
    by_exp = {}
    for r in phot:
        if expo[r["exposure"]]["object"] == "target_field" and r["source"] in cat:
            by_exp.setdefault(r["exposure"], []).append(r)
    for ex, rows in by_exp.items():
        det = np.array([[float(r["x_px"]) - x0, float(r["y_px"]) - y0] for r in rows])
        sky = np.array([cat[r["source"]] for r in rows])          # (east, north)
        # sky angle = detector angle + PA, with sky angle measured from north through east
        zs = sky[:, 1] + 1j * sky[:, 0]
        zd = det[:, 0] + 1j * det[:, 1]
        angle[ex] = float(np.degrees(np.angle(np.sum(zs * np.conj(zd)))))
    out = []
    for r in phot:
        e = expo[r["exposure"]]
        use = clean[(e["object"], r["source"])]
        if len(use) >= 3:
            o = growth_total(use, [float(r[f"o_r{k}"]) for k in use])
            eb = growth_total(use, [float(r[f"e_r{k}"]) for k in use])
        else:
            o, eb = float(r[f"o_r{use[-1]}"]), float(r[f"e_r{use[-1]}"])
        dx, dy = (float(r["x_px"]) - x0) / 1000.0, (float(r["y_px"]) - y0) / 1000.0
        pa = angle.get(r["exposure"], float(e["instrument_pa_deg"]))
        out.append(dict(obj=e["object"], src=r["source"], lr=np.log(o / eb), w=np.sqrt(o * eb / (o + eb)),
                        a=np.radians(4 * (float(e["hwp_deg"]) - phi0)), pa=np.radians(2 * pa),
                        rad=np.hypot(dx, dy) ** 2, ph=2 * np.arctan2(dy, dx)))
    return sheet, out


def fit_stokes(rows, eps, inst_extra=(0.0, 0.0), frame="sky"):
    """Fits (q, u) [sky or instrument frame] and one log transmission per source."""
    names = sorted({r["src"] for r in rows})
    idx = {n: i for i, n in enumerate(names)}

    def resid(p):
        q, u = p[0], p[1]
        t = p[4:] if frame == "both" else p[2:]
        out = []
        for r in rows:
            if frame == "sky":
                Q = q * np.cos(r["pa"]) + u * np.sin(r["pa"])
                U = -q * np.sin(r["pa"]) + u * np.cos(r["pa"])
            elif frame == "both":       # instrument-fixed (q, u) plus a sky-fixed (p[2], p[3]) turning with the rotator
                Q = q + p[2] * np.cos(r["pa"]) + p[3] * np.sin(r["pa"])
                U = u - p[2] * np.sin(r["pa"]) + p[3] * np.cos(r["pa"])
            else:
                Q, U = q, u
            m = eps * (Q + inst_extra[0]) * np.cos(r["a"]) + eps * (U + inst_extra[1]) * np.sin(r["a"])
            out.append((t[idx[r["src"]]] + np.log((1 + m) / (1 - m)) - r["lr"]) * r["w"])
        return np.array(out)

    sol = least_squares(resid, [0.0, 0.0] + ([0.0, 0.0] if frame == "both" else []) + [0.0] * len(names)).x
    return sol[0], sol[1]


def reduce(folder):
    sheet, rows = records(folder)
    sel = lambda obj, src=None, skip=None: [r for r in rows if r["obj"] == obj and (src is None or r["src"] == src) and r["src"] != skip]
    un = sel("unpolarised_standard")
    two = np.ptp([r["pa"] for r in un]) > 0.1      # two rotator settings separate the instrument from the star
    ip = fit_stokes(un, 1.0, frame="both" if two else "instrument")          # apparent (times eps)
    lam = sheet["filter_wavelength_um"][sheet["filter"]]
    cat = sheet["polarised_standard"]
    p_std = serkowski(cat["p_max_percent"] / 100.0, cat["lambda_max_um"], cat["k"], lam)
    pol = sel("polarised_standard")
    names = sorted({r["src"] for r in pol})

    def std_resid(p):
        eps, th, t = p[0], np.radians(p[1]), p[2]
        q, u = p_std * np.cos(2 * th), p_std * np.sin(2 * th)
        out = []
        for r in pol:
            Q = q * np.cos(r["pa"]) + u * np.sin(r["pa"])
            U = -q * np.sin(r["pa"]) + u * np.cos(r["pa"])
            m = eps * Q * np.cos(r["a"]) + eps * U * np.sin(r["a"]) + ip[0] * np.cos(r["a"]) + ip[1] * np.sin(r["a"])
            out.append((t + np.log((1 + m) / (1 - m)) - r["lr"]) * r["w"])
        return np.array(out)

    eps = least_squares(std_resid, [1.0, cat["angle_deg"], 0.0]).x[0]
    ip_true = (ip[0] / eps, ip[1] / eps)
    tgt = sheet["target"]
    field = sel("target_field", skip=tgt)
    # stage 1: each field star's instrument-frame Stokes (instrumental + its own sky polarisation rotated)
    stars = sorted({r["src"] for r in field})
    rowsA, vals = [], []
    for s in stars:
        rr = [r for r in field if r["src"] == s]
        qi, ui = fit_stokes(rr, eps, inst_extra=ip_true, frame="instrument")
        pa = np.mean([r["pa"] for r in rr])
        c, si = np.cos(pa), np.sin(pa)
        rad, ph = rr[0]["rad"], rr[0]["ph"]
        # stage 2: instrument-frame value = R(-PA) * common sky vector + k * rad * (cos ph, sin ph)
        rowsA += [[c, si, rad * np.cos(ph)], [-si, c, rad * np.sin(ph)]]
        vals += [qi, ui]
    k = np.linalg.lstsq(np.array(rowsA), np.array(vals), rcond=None)[0][2]
    t_rows = sel("target_field", tgt)
    rad, ph = t_rows[0]["rad"], t_rows[0]["ph"]
    extra = (ip_true[0] + k * rad * np.cos(ph), ip_true[1] + k * rad * np.sin(ph))
    q, u = fit_stokes(t_rows, eps, inst_extra=extra, frame="sky")
    return {"target_p_percent": float(100 * np.hypot(q, u)), "target_angle_deg": float((0.5 * np.degrees(np.arctan2(u, q))) % 180.0),
            "efficiency": float(eps)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(Path(sys.argv[1])), indent=1) + "\n")
