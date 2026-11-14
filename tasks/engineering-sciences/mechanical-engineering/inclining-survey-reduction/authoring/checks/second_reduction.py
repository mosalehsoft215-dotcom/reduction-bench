"""A second reduction of the inclining records, written without the reference program.

Differences from solution/incline_reduction.py: the hull is handled station by station (each offset section is a
polygon clipped by the inclined waterline, then Simpson's rule along the ship, with no interpolation between
stations); spaces open to the sea are taken out analytically; and the unknown centre of gravity is found by
simulating every condition (heel found by root search with the displaced volume held) and fitting the
simulated heels to the pendulum heels, instead of writing one moment balance per condition.

usage: python3 second_reduction.py <vessel folder>   -> prints a JSON line
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq, least_squares


def table(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def clip(poly, k, off):
    """keep the part of polygon [(y, z)] with z - k*y <= off (Sutherland-Hodgman)"""
    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        fp, fq = p[1] - k * p[0] - off, q[1] - k * q[0] - off
        if fp <= 0:
            out.append(p)
        if fp * fq < 0:
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def area_props(poly):
    if len(poly) < 3:
        return 0.0, 0.0, 0.0
    y = np.array([p[0] for p in poly])
    z = np.array([p[1] for p in poly])
    y2, z2 = np.roll(y, -1), np.roll(z, -1)
    cr = y * z2 - y2 * z
    a = cr.sum() / 2
    if abs(a) < 1e-12:
        return 0.0, 0.0, 0.0
    return a, ((y + y2) * cr).sum() / (6 * a), ((z + z2) * cr).sum() / (6 * a)


class Stations:
    def __init__(self, folder, part):
        off = table(folder / "offsets.csv")
        self.xs = np.array(sorted({float(r["x_m"]) for r in off}))
        self.zs = np.array(sorted({float(r["z_m"]) for r in off}))
        hb = {(float(r["x_m"]), float(r["z_m"])): float(r["half_breadth_m"]) for r in off}
        self.hb = np.array([[hb[(x, z)] for z in self.zs] for x in self.xs])
        self.polys = []
        for i in range(self.xs.size):
            b = self.hb[i]
            stb = [(b[j], self.zs[j]) for j in range(self.zs.size) if b[j] > 0 or j == self.zs.size - 1]
            self.polys.append(stb + [(-y, z) for y, z in reversed(stb)])
        self.spaces = part.get("open_to_sea", [])

    def width(self, x, z):
        col = np.array([np.interp(z, self.zs, self.hb[i]) for i in range(self.xs.size)])
        return float(np.interp(x, self.xs, col))

    def immersed(self, a, bx, cy):
        A, Y, Z = [], [], []
        for x, poly in zip(self.xs, self.polys):
            ar, yc, zc = area_props(clip(poly, cy, a + bx * x)) if len(poly) >= 3 else (0, 0, 0)
            A.append(ar)
            Y.append(ar * yc)
            Z.append(ar * zc)
        A, Y, Z = map(np.array, (A, Y, Z))
        v = simpson(A, x=self.xs)
        m = np.array([simpson(A * self.xs, x=self.xs), simpson(Y, x=self.xs), simpson(Z, x=self.xs)])
        for s in self.spaces:
            if s["kind"] == "moonpool":     # prism from the bottom up to the inclined waterline
                L, Bm = s["length_m"], s["breadth_m"]
                h = a + bx * s["x_m"] + cy * s["y_m"]
                vol = L * Bm * h
                mx = vol * s["x_m"] + bx * L ** 3 * Bm / 12
                my = vol * s["y_m"] + cy * Bm ** 3 * L / 12
                mz = L * Bm * h * h / 2 + (bx ** 2 * L ** 3 * Bm + cy ** 2 * Bm ** 3 * L) / 24
            else:                            # tunnel wholly submerged, side to side at its axis height
                r = s["diameter_m"] / 2
                vol = np.pi * r * r * 2 * self.width(s["x_m"], s["z_m"])
                mx, my, mz = vol * s["x_m"], 0.0, vol * s["z_m"]
            v -= vol
            m -= np.array([mx, my, mz])
        return v, m


def box_liquid(t, a, bx, cy, n=240):
    gx = t["xa"] + (np.arange(n) + 0.5) * (t["xf"] - t["xa"]) / n
    gy = t["yp"] + (np.arange(n) + 0.5) * (t["ys"] - t["yp"]) / n
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    h = np.clip(a + bx * X + cy * Y - t["zb"], 0, t["zt"] - t["zb"])
    dA = (t["xf"] - t["xa"]) * (t["ys"] - t["yp"]) / n / n
    return h.sum() * dA, np.array([(h * X).sum(), (h * Y).sum(), (h * (t["zb"] + h / 2)).sum()]) * dA


def reduce(folder):
    folder = Path(folder)
    part = json.loads((folder / "ship.json").read_text())
    hull = Stations(folder, part)
    rho = np.median([float(r["density_kg_m3"]) for r in table(folder / "water_samples.csv")]) / 1000

    marks = {m["mark"]: m["x_m"] for m in part["draft_marks"]}
    reads = {}
    for r in table(folder / "draft_readings.csv"):
        reads.setdefault((r["mark"], r["side"]), []).append(float(r["reading_m"]))
    pts = []
    for (mk, side), v in reads.items():
        z = np.mean(v)
        pts.append((marks[mk], hull.width(marks[mk], z) * (1 if side == "starboard" else -1), z))
    pts = np.array(pts)
    # leave each mark scale out in turn; the one the others predict worst is dropped while it misses by > 6 cm
    while len(pts) > 4:
        miss = []
        for i in range(len(pts)):
            o = np.arange(len(pts)) != i
            cf = np.linalg.lstsq(np.c_[np.ones(o.sum()), pts[o, :2]], pts[o, 2], rcond=None)[0]
            miss.append(abs(pts[i, 2] - (cf[0] + cf[1] * pts[i, 0] + cf[2] * pts[i, 1])))
        worst = int(np.argmax(miss))
        if miss[worst] < 0.06:
            break
        pts = np.delete(pts, worst, axis=0)
    a0, bx, cy0 = np.linalg.lstsq(np.c_[np.ones(len(pts)), pts[:, :2]], pts[:, 2], rcond=None)[0]

    app_v = sum(a["volume_m3"] for a in part["appendages"])
    app_m = sum(a["volume_m3"] * np.array([a["x_m"], a["y_m"], a["z_m"]]) for a in part["appendages"])
    v0, _ = hull.immersed(a0, bx, cy0)
    disp = rho * (v0 + app_v)

    tanks = {}
    for r in table(folder / "tank_soundings.csv"):
        t = dict(xa=float(r["x_aft_m"]), xf=float(r["x_fwd_m"]), yp=float(r["y_port_m"]), ys=float(r["y_starboard_m"]),
                 zb=float(r["z_bottom_m"]), zt=float(r["z_top_m"]), rho=float(r["density_t_m3"]), remark=r["remark"])
        s = float(r["sounding_m"])
        if s <= 0:
            continue
        if "pressed full" in t["remark"]:
            t["vol"] = (t["xf"] - t["xa"]) * (t["ys"] - t["yp"]) * (t["zt"] - t["zb"])
        else:
            lev = t["zb"] + s - bx * float(r["pipe_x_m"]) - cy0 * float(r["pipe_y_m"])
            t["vol"] = box_liquid(t, lev, bx, cy0)[0]
        tanks[r["tank"]] = t
    sets, used = [], set()
    for name, t in tanks.items():
        if name in used:
            continue
        mates = [name]
        if t["remark"].startswith("cross-over valve to ") and t["remark"].endswith(" open"):
            other = t["remark"][len("cross-over valve to "):-len(" open")]
            if other in tanks:
                mates.append(other)
        used.update(mates)
        sets.append([tanks[m] for m in mates])

    def liquid(b, c):
        mass, mom = 0.0, np.zeros(3)
        for grp in sets:
            want = sum(t["vol"] for t in grp)
            if all("pressed full" in t["remark"] for t in grp):
                lev = 1e4
            else:
                f = lambda lev: sum(box_liquid(t, lev, b, c, 120)[0] for t in grp) - want
                lev = brentq(f, -60, 80, xtol=1e-7)
            for t in grp:
                v, m = box_liquid(t, lev, b, c)
                mass += t["rho"] * v
                mom += t["rho"] * m
        return mass, mom

    weights = {r["weight"]: r for r in table(folder / "test_weights.csv")}
    for r in table(folder / "load_cell.csv"):
        weights[r["weight"]]["mass_t"] = r["load_cell_t"]
    out_items, fit_items = [], []
    for r in table(folder / "lightship_list.csv"):
        m = float(r["mass_t"])
        if r["status"] == "fit":
            fit_items.append((m, np.array([r["x_m"], r["y_m"], r["z_m"]], float)))
        elif r["hung_from_z_m"]:
            out_items.append((m, np.array([r["hung_from_x_m"], r["hung_from_y_m"], r["hung_from_z_m"]], float)))
        else:
            out_items.append((m, np.array([r["x_m"], r["y_m"], r["z_m"]], float)))
    shifts = table(folder / "arrangements.csv")

    pend = table(folder / "battens.csv")
    tclk = np.array([float(r["clock_s"]) for r in pend])
    names = [p["pendulum"] for p in part["pendulums"]]
    lengths = {p["pendulum"]: p["length_m"] * 1000 for p in part["pendulums"]}
    steps = {}
    for n in names:
        v = []
        for s in shifts:
            sel = (tclk >= float(s["from_clock_s"])) & (tclk <= float(s["to_clock_s"]))
            v.append(np.median([float(pend[i][n]) for i in np.flatnonzero(sel)]) / lengths[n])
        steps[n] = np.array(v) - v[0]
    inc = table(folder / "inclinometer_log.csv")
    ic = np.array([float(r["clock_s"]) for r in inc])
    ia = np.array([float(r["heel_deg"]) for r in inc])
    v = [np.tan(np.radians(np.median(ia[(ic >= float(s["from_clock_s"])) & (ic <= float(s["to_clock_s"]))]))) for s in shifts]
    steps["inclinometer"] = np.array(v) - v[0]
    # a pendulum whose heel steps differ from the median of all three instruments by more than 3 % is dropped
    med = np.median(np.array(list(steps.values())), axis=0)
    big = np.abs(med) > 0.3 * np.abs(med).max()
    good = [n for n in names if np.median(steps[n][big] / med[big]) - 1 < 0.03 and 1 - np.median(steps[n][big] / med[big]) < 0.03]
    heel = np.mean([steps[n] for n in (good or ["inclinometer"])], axis=0)
    # arrangement 0 is repeated during the day; its change in time is drift, removed as a straight line
    tm = np.array([0.5 * (float(s["from_clock_s"]) + float(s["to_clock_s"])) for s in shifts])
    rep = [k for k, s in enumerate(shifts) if all(s[w] == shifts[0][w] for w in "ABCD")]
    slope = np.polyfit(tm[rep], heel[rep], 1)[0] if len(rep) > 1 else 0.0
    heel = cy0 + heel - slope * (tm - tm[0])

    def known(k, b, c):
        mass, mom = liquid(b, c)
        for m, p in out_items:
            mass += m
            mom = mom + m * p
        for w in weights.values():
            y = float(w["y_port_m"]) if shifts[k][w["weight"]] == "P" else float(w["y_starboard_m"])
            mass += float(w["mass_t"])
            mom = mom + float(w["mass_t"]) * np.array([float(w["x_m"]), y, float(w["z_m"])])
        return mass, mom

    cache = {}

    def level_for(c):
        key = round(c, 9)
        if key not in cache:
            g = lambda a: hull.immersed(a, bx, c)[0] - v0
            a = brentq(g, a0 - 0.5, a0 + 0.5, xtol=1e-8)
            v, m = hull.immersed(a, bx, c)
            cache[key] = (v + app_v, (m + app_m) / (v + app_v))
        return cache[key]

    def balance(c, k, yz):
        vb, cb = level_for(c)
        mass, mom = known(k, bx, c)
        mf = disp - mass
        g = (mom + mf * np.array([0.0, yz[0], yz[1]])) / disp
        return (g[1] + c * g[2]) - (cb[1] + c * cb[2])

    def simulate(yz):
        out = []
        for k in range(len(shifts)):
            out.append(brentq(lambda c: balance(c, k, yz), heel[k] - 0.05, heel[k] + 0.05, xtol=1e-10))
        return np.array(out)

    mass0, mom0 = known(0, bx, cy0)
    mf = disp - mass0
    _, cb0 = level_for(cy0)
    # starting point: the balance at the measured heels is linear in (y, z)
    A, r = [], []
    for k in range(len(shifts)):
        b0 = balance(heel[k], k, (0.0, 0.0))
        A.append([balance(heel[k], k, (1.0, 0.0)) - b0, balance(heel[k], k, (0.0, 1.0)) - b0])
        r.append(-b0)
    start = np.linalg.lstsq(np.array(A), np.array(r), rcond=None)[0]
    fit = least_squares(lambda yz: (simulate(yz) - heel) * 1e4, x0=start, x_scale=[0.1, 0.5])
    yf, zf = fit.x
    xf = (disp * (cb0[0] + bx * cb0[2]) - (mom0[0] + bx * mom0[2])) / mf - bx * zf
    m_ls = mf + sum(m for m, _ in fit_items)
    g = (mf * np.array([xf, yf, zf]) + sum(m * p for m, p in fit_items)) / m_ls
    return dict(light_mass_t=float(m_ls), light_lcg_m=float(g[0]), light_vcg_m=float(g[2]))


if __name__ == "__main__":
    print(json.dumps(reduce(sys.argv[1])))
