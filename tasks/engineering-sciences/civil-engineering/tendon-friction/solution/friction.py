"""Duct friction of a post-tensioned girder (mu per radian, K per metre) and the force of one tendon at one station,
from one bridge's stressing records.

    python3 friction.py <bridge folder> <output json>

Route.
1. Every duct is rebuilt as a curve in space: the girder axis in plan (straights and arcs by station), the soffit
   level from the vertical alignment, the height of the duct centre above the soffit (one parabola in station between
   two control points, measured vertically) and its lateral place: a web duct in the mid-plane of its web, which leans
   outward, so its offset from the axis grows with its height; a bottom-slab duct at its stated offset. The strands
   do not follow the duct centre: tensioned, they take the shortest line the duct allows, within the bundle's
   eccentricity of the centre (centred at both anchor heads), bearing on the inside of each bend and crossing the
   duct on a straight chord where a bend reverses or a straight runs into a bend. That line is solved on a 0.1 m
   station grid; the friction angle is the turn of its tangent from anchorage A, the distance its arc length.
2. Gauge readings become jack forces through each jack's own calibration table (quadratic through the table, in the
   gauge's own unit). The dead-end load cells, less their stage-0 zero readings, give ln(P_dead / P_jack) =
   -(mu theta + K s) per tendon and stage; mu and K are their least-squares fit. A bridge without load cells gets
   mu and K from the ram travels alone (step 3's prediction fitted to every jacking's travels, each jacking with its
   own offset, since the 20 % mark is not where the strands start to move).
3. Before a cell is trusted, every tendon's ram travel is predicted with the same mu and K (stage forces from the
   gauges, travel counted from the 20 % mark, strand inside the jack, the force already in a tendon stressed from its
   other end, the draw-in at each lock-off: the marks sit behind the jack's gripping wedges, so their movement is
   the draw-in plus the recoil of the free strand between the anchor head and the grip, F (L + e) / EA, where L is
   the length with the ram closed and e the ram's extension since its last grip (the strand drawn out of the duct
   since then sits in the jack too), which is taken off; and the restrokes: the
   ram starts closed, and at each stage hold where the next stage would not fit in the stroke left the anchor wedges
   seat by that end's draw-in and the strand in the jack goes slack, so the next stroke first travels the draw-in
   plus F L / EA before the force rises again, and that travel is in the added-up readings): a load-cell tendon whose
   travel falls short of the others', relative to the prediction, has a duct that loses more than mu and K say (a
   local kink), and its cell is left out;
   so is a cell whose tendon's travel falls short while the remaining cells agree with one mu and K to their own
   precision and all of them together clearly do not.
4. The asked-for force follows from the same model of that tendon's stressing and lock-offs. If that tendon is itself
   one with a damaged duct, its local loss is part of the force it carries: the size of the loss comes from its own
   load cell (the excess of -ln(P_dead / P_jack) over mu theta + K s), its station from its own ram travel.
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import splu

DS = 0.02
OPTIONS = dict(zero=True, outliers=True, plan="magnitude", strand_in_jack=True, mark=True, own_jack=True, units=True,
               lean=True, vertical=True, own_loss=True, bundle=True, recoil=True, ram_length=True,
               nominal=False, restroke=True)


def read(folder):
    folder = Path(folder)
    b = json.loads((folder / "bridge.json").read_text())
    knots = {}
    for r in csv.DictReader(open(folder / "tendon_profiles.csv")):
        knots.setdefault(r["tendon"], []).append((float(r["station_m"]), float(r["height_m"]), float(r["slope"])))
    rec = {}
    for r in csv.DictReader(open(folder / "stressing_records.csv")):
        rec.setdefault(r["tendon"], []).append(r)
    return b, knots, rec


def cumtrapz(f, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def height(knots, st):
    """height above the soffit: one parabola between two control points (heights and slopes at both ends)"""
    k = sorted(knots)
    h = np.zeros_like(st)
    for (x0, h0, s0), (x1, _, s1) in zip(k[:-1], k[1:]):
        m = (st >= x0) & (st <= x1)
        u = st[m] - x0
        h[m] = h0 + s0 * u + 0.5 * (s1 - s0) / (x1 - x0) * u * u
    return h


def soffit(b, st, o):
    v = b["girder"]["vertical_alignment"]
    if not isinstance(v, dict) or not o["vertical"]:
        return np.zeros_like(st)
    a, e, r = v["from_station_m"], v["to_station_m"], v["radius_m"]
    sign = -1.0 if v["geometry"].startswith("crest") else 1.0
    u = np.clip(st, a, e) - a
    return sign * (u * u / (2.0 * r) + np.where(st > e, (st - e) * (e - a) / r, 0.0))


def lateral(b, tendon, h, o):
    pos = next(t for t in b["tendons"] if t["name"] == tendon)["position"]
    sec = b["girder"]["cross_section"]
    if "web" in pos and not pos.startswith("bottom"):
        side = 1.0 if pos.startswith("left") else -1.0
        lean = sec["web_lean_outward_horizontal_per_vertical"] if o["lean"] else 0.0
        if not o["lean"]:
            h = np.full_like(h, np.mean(h))
        return side * (sec["web_mid_plane_offset_at_soffit_m"] + h * lean)
    m = re.search(r"([\d.]+) m (left|right)", pos)
    return np.full_like(h, float(m.group(1)) * (1.0 if m.group(2) == "left" else -1.0))


def angle(b, tendon, knots, o):
    """station grid, arc length and the angle turned through from A"""
    length = b["girder"]["length_between_anchorages_m"]
    st = np.arange(0.0, length + DS / 2, DS)
    kp = np.zeros_like(st)
    for seg in b["girder"]["plan_alignment"]:
        if seg["geometry"] == "circular arc":
            sign = 1.0 if seg["turn"] == "left" else -1.0
            if o["plan"] == "none":
                sign = 0.0
            kp[(st >= seg["from_station_m"]) & (st < seg["to_station_m"])] = sign / seg["radius_m"]
    psi = cumtrapz(kp, st)
    px, py = cumtrapz(np.cos(psi), st), cumtrapz(np.sin(psi), st)
    h = height(knots, st)
    y = lateral(b, tendon, h, o)
    pts = np.stack([px - y * np.sin(psi), py + y * np.cos(psi), soffit(b, st, o) + h], 1)
    ecc = b.get("ducts", {}).get("strand_bundle_eccentricity_mm", 0.0) / 1e3
    if not o["bundle"] or ecc <= 0.0:
        return (st,) + turned(pts)
    sub = np.arange(0, st.size, 5)
    if sub[-1] != st.size - 1:
        sub = np.append(sub, st.size - 1)
    s_sub, th_sub = turned(bundle_path(pts[sub], ecc))
    return st, np.interp(st, st[sub], s_sub), np.interp(st, st[sub], th_sub)


def turned(pts):
    """arc length and turned angle of a polyline, the turn between two chords placed at their common node"""
    d = np.diff(pts, axis=0)
    seg = np.linalg.norm(d, axis=1)
    t = d / seg[:, None]
    turn = np.arctan2(np.linalg.norm(np.cross(t[1:], t[:-1]), axis=1), np.sum(t[1:] * t[:-1], 1))
    return np.concatenate([[0.0], np.cumsum(seg)]), np.concatenate([[0.0, 0.0], np.cumsum(turn)])


def bundle_path(c, ecc, rho=0.05, tol=1e-10):
    """where the tensioned strands run: the shortest line through the duct, i.e. inside the tube of radius ecc
    around the duct centre c, centred at both anchor heads. Each node may move within the plane normal to the duct;
    the sum of squared chords is minimised under the disc constraint by ADMM (sparse normal equations, projection)"""
    n = len(c)
    t = np.gradient(c, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    e1 = np.cross([0.0, 0.0, 1.0], t)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(t, e1)
    i = np.arange(n - 1)
    rows, cols, vals = [], [], []
    for k in range(3):
        for node, sign in ((i + 1, 1.0), (i, -1.0)):           # chord i = c[i+1] - c[i] + offset[i+1] - offset[i]
            ok = (node >= 1) & (node <= n - 2)
            for col, e in ((0, e1), (1, e2)):
                rows.append(3 * i[ok] + k)
                cols.append(2 * (node[ok] - 1) + col)
                vals.append(sign * e[node[ok], k])
    a = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(3 * (n - 1), 2 * (n - 2)))
    rhs0 = -(a.T @ np.diff(c, axis=0).ravel())
    solve = splu((a.T @ a + rho * sparse.identity(2 * (n - 2))).tocsc()).solve
    z, w = np.zeros(2 * (n - 2)), np.zeros(2 * (n - 2))
    for it in range(20000):
        u = solve(rhs0 + rho * (z - w))
        v = (u + w).reshape(-1, 2)
        z_new = (v * np.minimum(1.0, ecc / np.maximum(np.linalg.norm(v, axis=1), 1e-15))[:, None]).ravel()
        w += u - z_new
        settled = it > 50 and np.abs(z_new - z).max() < tol and np.abs(u - z_new).max() < tol
        z = z_new
        if settled:
            break
    off = np.zeros((n, 2))
    off[1:-1] = z.reshape(-1, 2)
    return c + off[:, :1] * e1 + off[:, 1:] * e2


def calibration(jack, o):
    g, f = np.array(jack["calibration"]["gauge"], float), np.array(jack["calibration"]["force_kN"], float)
    if not o["units"] and jack["gauge_unit"] == "psi":
        g = g * 0.0689476                      # table read as if the gauge were in bar
    c = np.polyfit(g, f, 2)
    return lambda p: np.polyval(c, p) * 1e3


def campaign(b, rec, o, ea):
    """per tendon and end: stage forces (N), relative ram travel (m), load-cell forces (N), draw-in (m)"""
    jacks = {end: b["jacks"][end if o["own_jack"] else "A"] for end in ("A", "B")}
    cal = {end: calibration(j, o) for end, j in jacks.items()}
    out = {}
    for t, rows in rec.items():
        ends = []
        for end in dict.fromkeys(r["end"] for r in rows):
            rr = [r for r in rows if r["end"] == end]
            st = [r for r in rr if r["stage"] not in ("0", "lock-off")]
            zero = [float(r["dead_end_load_cell_kN"]) for r in rr if r["stage"] == "0"]
            force = np.array([cal[end](float(r["gauge"])) for r in st])
            if o["nominal"]:
                force = np.array([float(r["stage"]) / 100.0 * b["jacking_force_kN"] * 1e3 for r in st])
            travel = np.array([float(r["extension_mm"]) for r in st]) / 1e3
            cell = np.array([float(r["dead_end_load_cell_kN"]) if r["dead_end_load_cell_kN"] else np.nan for r in st]) * 1e3
            if zero and o["zero"]:
                cell = cell - zero[0] * 1e3
            draw = [float(r["draw_in_mm"]) / 1e3 for r in rr if r["stage"] == "lock-off"][0]
            if o["recoil"]:
                # the marks sit behind the jack's gripping wedges: when the jack lets go, the strand inside it
                # shortens elastically by F L / EA, and the marks move by that as well as by the draw-in
                draw -= force[-1] * jacks[end]["strand_length_in_jack_m"] / ea
            ends.append(dict(end=end, force=force, travel=travel, cell=cell, draw=draw,
                             in_jack=jacks[end]["strand_length_in_jack_m"] if o["strand_in_jack"] else 0.0,
                             stroke=jacks[end].get("stroke_mm", 1e9) / 1e3 if o["restroke"] else 1e9))
        out[t] = ends
    return out


def jack_curve(force, s, theta, mu, k, end):
    if end == "A":
        return force * np.exp(-(mu * theta + k * s))
    return force * np.exp(-(mu * (theta[-1] - theta) + k * (s[-1] - s)))


def lock_off(p, draw, ea, s, end):
    q = p if end == "A" else p[::-1]
    u = s if end == "A" else (s[-1] - s)[::-1]
    slack = cumtrapz(q, u) - q ** 2 * cumtrapz(1.0 / q, u)
    hit = np.nonzero(slack[1:] >= draw * ea)[0]
    i = hit[0] + 1 if hit.size else q.size - 1
    out = q.copy()
    out[: i + 1] = q[i] ** 2 / q[: i + 1]
    return out if end == "A" else out[::-1]


def replay(ends, s, theta, mu, k, ea, o):
    """predicted relative ram travel per end and stage, and the force profile after the last lock-off.
    The ram starts closed at zero force and is not run out of stroke between stages: at a stage hold, if the
    travel to the next stage would not fit in what is left of the stroke, the anchor wedges take the load at that
    stage's force F and draw in, the strand inside the jack is unloaded, and after the new grip the ram first travels
    the draw-in plus F L / EA before the strands move on; that travel is part of the added-up readings."""
    existing, travel = np.zeros_like(s), []
    for e in ends:
        start = existing

        def total(f):
            new = np.maximum(start, jack_curve(f, s, theta, mu, k, e["end"]))
            drawn = float(np.sum(0.5 * ((new - start)[1:] + (new - start)[:-1]) * np.diff(s)))
            return drawn / ea + f * e["in_jack"] / ea

        pred = np.array([total(f) for f in e["force"]])
        draw = e["draw"]
        for _ in range(3 if o["ram_length"] and o["recoil"] else 1):
            extra, ram, regripped = np.zeros_like(pred), 0.0, False
            for i in range(len(pred)):
                ram += pred[i] - (pred[i - 1] if i else 0.0)
                if i + 1 < len(pred) and ram + pred[i + 1] - pred[i] > e["stroke"]:
                    ram = max(draw, 0.0) + e["force"][i] * e["in_jack"] / ea
                    extra[i + 1:] += ram
                    regripped = True
            if o["ram_length"] and o["recoil"]:
                # at lock-off the free strand between the anchor head and the gripping wedges is the closed-ram
                # length plus the ram's extension since its last grip (which holds the draw-in after a restroke):
                # marks = draw-in + F (L + ram) / EA, and e["draw"] already has F L / EA taken off
                eps = e["force"][-1] / ea
                rest = ram - (max(draw, 0.0) if regripped else 0.0)
                draw = (e["draw"] - eps * rest) / (1.0 + (eps if regripped else 0.0))
        if o["mark"]:
            pred = pred - pred[0] + extra - extra[0]
        else:
            pred = pred + extra
        travel.append(pred)
        existing = np.maximum(existing, jack_curve(e["force"][-1], s, theta, mu, k, e["end"]))
        existing = lock_off(existing, draw, ea, s, e["end"])
    return travel, existing


def fit_travels(data, geo, ea, o, start=(0.2, 1.0e-3)):
    """mu and K from the ram travels of every jacking alone (bridges without load cells)"""
    from scipy.optimize import least_squares

    def resid(v):
        out = []
        for t, ends in data.items():
            pred, _ = replay(ends, *geo[t][1:], v[0], v[1] * 1e-3, ea, o)
            for e, p in zip(ends, pred):
                # each jacking gets its own offset: the force at the 20 % mark is read with a gauge error
                r = e["travel"] - p
                out.extend((r - r.mean()) * 1e3)
        return np.array(out)

    v = least_squares(resid, [start[0], start[1] * 1e3], x_scale=[0.02, 0.1], diff_step=1e-4).x
    return v[0], v[1] * 1e-3


def fit_cells(cells):
    """least squares of -ln(P_dead / P_jack) = mu theta + K s over every cell reading"""
    a = np.array([[th, ln] for th, ln, _ in cells])
    y = np.array([v for _, _, v in cells])
    return np.linalg.lstsq(a, y, rcond=None)[0]


def own_loss(ends, st, s, theta, mu, k, ea, o):
    """extra angle and station of a damaged duct from its own load cell (size) and ram travel (place)"""
    e = ends[0]
    ok = np.isfinite(e["cell"]) & (e["force"] > 0.3 * e["force"].max())
    extra = np.mean(-np.log(e["cell"][ok] / e["force"][ok]) - (mu * theta[-1] + k * s[-1])) / mu
    best = None
    for where in np.arange(1.0, st[-1] - 1.0, 0.5):
        th = theta + extra * (st >= where)
        pred, _ = replay(ends[:1], s, th, mu, k, ea, o)
        err = float(np.sum((pred[0][1:] - e["travel"][1:]) ** 2))
        if best is None or err < best[0]:
            best = (err, where)
    return extra, best[1]


def solve(folder, **change):
    o = dict(OPTIONS, **change)
    b, knots, rec = read(folder)
    lot = b["strand_lot"]
    ea = lot["strands_per_tendon"] * lot["strand_area_mm2"] * 1e-6 * lot["modulus_GPa"] * 1e9
    data = campaign(b, rec, o, ea)
    geo = {t: angle(b, t, knots[t], o) for t in data}
    cells = {}
    for t, ends in data.items():
        e = ends[0]
        if len(ends) == 1 and np.isfinite(e["cell"]).any():
            _, s, th = geo[t]
            ok = np.isfinite(e["cell"]) & (e["force"] > 0.3 * e["force"].max())
            cells[t] = [(th[-1], s[-1], -np.log(c / f)) for c, f in zip(e["cell"][ok], e["force"][ok])]

    def deviations(mu, k):
        """relative shortfall of the full-stage ram travel at each tendon's first jacking end"""
        out = {}
        for t, ends in data.items():
            pred, _ = replay(ends[:1], *geo[t][1:], mu, k, ea, o)
            out[t] = (ends[0]["travel"][-1] - pred[0][-1]) / pred[0][-1]
        return out

    def cell_rms(group, mu, k):
        return float(np.sqrt(np.mean([(v - (mu * th + k * ln)) ** 2 for t in group for th, ln, v in cells[t]])))

    used = set(cells)
    if used:
        mu, k = fit_cells([c for t in used for c in cells[t]])
    else:
        mu, k = fit_travels(data, geo, ea, o)
    dev = deviations(mu, k)
    while o["outliers"] and len(used) > 2:
        # leave each load-cell tendon out in turn: a damaged duct shows as a travel shortfall far outside the scatter
        # of all other tendons once mu and K no longer contain its extra loss, and the other cells then agree with
        # one mu and K to their own precision while all of them together do not
        now = cell_rms(used, mu, k)
        trials = []
        for c in used:
            m, kk = fit_cells([v for t in used - {c} for v in cells[t]])
            d = deviations(m, kk)
            rest = np.array([v for t, v in d.items() if t != c])
            mid, spread = np.median(rest), 1.4826 * np.median(np.abs(rest - np.median(rest)))
            trials.append(((d[c] - mid) / spread, d[c] - mid, c, m, kk, d, cell_rms(used - {c}, m, kk)))
        z, gap, c, m, kk, d, rms_rest = min(trials, key=lambda r: r[0])
        travel_says = z <= -8.0 and gap <= -0.008
        cells_say = gap <= -0.005 and now > 0.005 and now > 3.0 * rms_rest
        if not (travel_says or cells_say):
            break
        used, mu, k, dev = used - {c}, m, kk, d
    tendon, station = b["force_asked_for"]["tendon"], b["force_asked_for"]["station_m"]
    st, s, theta = geo[tendon]
    damage = None
    if o["own_loss"] and tendon in set(cells) - used:
        extra, where = own_loss(data[tendon], st, s, theta, mu, k, ea, o)
        theta = theta + extra * (st >= where)
        damage = dict(extra_angle_rad=round(float(extra), 4), station_m=float(where))
    _, final = replay(data[tendon], s, theta, mu, k, ea, o)
    force = float(np.interp(station, st, final)) / 1e3
    resid = [v - (mu * th + k * ln) for t in used for th, ln, v in cells[t]] or [0.0]
    info = dict(cell_rms_permille=round(1e3 * float(np.sqrt(np.mean(np.square(resid)))), 3), cells_used=sorted(used), cells_left_out=sorted(set(cells) - used), damaged_asked_tendon=damage,
                travel_deviation_percent={t: round(100 * v, 3) for t, v in dev.items()})
    return dict(mu=float(mu), K_per_m=float(k), force_kN=force), info


if __name__ == "__main__":
    result, info = solve(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(dict(result, details=info), indent=1) + "\n")
