#!/usr/bin/env python3
"""Reduces one inclining survey to the lightship mass and centre of gravity.

usage: python3 incline_reduction.py <survey folder> <result.json>

Route: body from the offsets -> waterplane from the draft scales -> displacement and centre of float_at ->
tank liquid_moments from the soundings -> heel of every condition from the pendulums -> the unknown part of the
ship from the equilibrium of all conditions (exact centres of float_at and liquid at each heel) ->
lightship by removing what does not belong and fitting what was missing.
"""
import csv
import json
import os
import re
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

SWITCHES = dict(pipe_plane=True, exact_liquid=True, liquid_moves=True, hang_point=True, open_spaces=True,
               exact_heel=True, trim_lever=True, mark_positions=True, joined_tanks=True, well_waterplane=True,
               load_cell=True, sample_screen=True, mark_screen=True, pendulum_screen=True, drift=True, inclinometer_only=False)


def read_table(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


class HullFromOffsets:
    def __init__(self, folder, part, o, nx=1600, ny=900, dz=0.01):
        off = read_table(os.path.join(folder, "offsets.csv"))
        xs = np.array(sorted({float(r["x_m"]) for r in off}))
        zs = np.array(sorted({float(r["z_m"]) for r in off}))
        hb = np.zeros((xs.size, zs.size))
        ix = {v: i for i, v in enumerate(xs)}
        iz = {v: j for j, v in enumerate(zs)}
        for r in off:
            hb[ix[float(r["x_m"])], iz[float(r["z_m"])]] = float(r["half_breadth_m"])
        ends = read_table(os.path.join(folder, "profile.csv"))
        ez = np.array([float(r["z_m"]) for r in ends])
        ea = np.array([float(r["aft_end_x_m"]) for r in ends])
        ef = np.array([float(r["fwd_end_x_m"]) for r in ends])
        self.D = float(part["depth_m"])
        self.B = float(part["beam_m"])
        x0, x1 = xs.min(), xs.max()
        self.dx = (x1 - x0) / nx
        self.xg = x0 + (np.arange(nx) + 0.5) * self.dx
        self.dy = self.B / ny
        self.yg = -self.B / 2 + (np.arange(ny) + 0.5) * self.dy
        z_rest = np.arange(0.0, self.D + dz / 2, dz)
        # each station in z, then each waterline in x through the body ends
        hz = np.array([PchipInterpolator(zs, hb[i])(z_rest) for i in range(xs.size)])
        aft = np.interp(z_rest, ez, ea)
        fwd = np.interp(z_rest, ez, ef)
        table = np.zeros((nx, z_rest.size))
        for j in range(z_rest.size):
            keep = (xs > aft[j] + 1e-6) & (xs < fwd[j] - 1e-6)
            px = np.concatenate([[aft[j]], xs[keep], [fwd[j]]])
            py = np.concatenate([[0.0], np.maximum(hz[keep, j], 0.0), [0.0]])
            v = PchipInterpolator(px, py)(np.clip(self.xg, aft[j], fwd[j]))
            table[:, j] = np.where((self.xg < aft[j]) | (self.xg > fwd[j]), 0.0, np.maximum(v, 0.0))
        table = np.maximum.accumulate(table, axis=1)          # sides do not fall in
        self.table, self.zf = table, z_rest
        zb = np.full((nx, ny), np.inf)
        ay = np.abs(self.yg)
        ramp = np.arange(z_rest.size) * 1e-9
        for i in range(nx):
            col = table[i] + ramp
            inside = ay <= col[-1]
            zb[i, inside] = np.interp(ay[inside], col, z_rest)
        self.X, self.Y = np.meshgrid(self.xg, self.yg, indexing="ij")
        self.mask = np.isfinite(zb)
        self.spaces = part.get("open_to_sea", []) if o["open_spaces"] else []
        self.wells = []
        for s in self.spaces:
            if s["kind"] == "moonpool" and not o["well_waterplane"]:
                self.wells.append(s)       # taken off as a fixed upright volume only (misreading study)
            elif s["kind"] == "moonpool":
                m = (np.abs(self.X - s["x_m"]) <= s["length_m"] / 2) & (np.abs(self.Y - s["y_m"]) <= s["breadth_m"] / 2)
                self.mask &= ~m
        self.zb = np.where(self.mask, zb, self.D)

    def half_breadth(self, x, z):
        i = np.clip(np.round((x - self.xg[0]) / self.dx).astype(int), 0, self.xg.size - 1)
        return np.interp(z, self.zf, self.table[i])

    def immersed(self, a, slope_x, slope_y):
        top = np.minimum(a + slope_x * self.X + slope_y * self.Y, self.D)
        h = np.where(self.mask, np.clip(top - self.zb, 0.0, None), 0.0)
        lo, hi = self.zb, self.zb + h
        hz = h * 0.5 * (lo + hi)
        for s in self.spaces:
            if s["kind"] == "tunnel":
                r = s["diameter_m"] / 2
                d = np.sqrt(np.maximum(r * r - (self.X - s["x_m"]) ** 2, 0.0))
                g0 = np.maximum(s["z_m"] - d, lo)
                g1 = np.minimum(s["z_m"] + d, hi)
                g = np.where(self.mask & (d > 0), np.clip(g1 - g0, 0.0, None), 0.0)
                h = h - g
                hz = hz - g * 0.5 * (g0 + g1)
        dA = self.dx * self.dy
        v = h.sum() * dA
        m = np.array([(h * self.X).sum(), (h * self.Y).sum(), hz.sum()]) * dA
        for s in self.wells:
            d = a + slope_x * s["x_m"]
            w = s["length_m"] * s["breadth_m"] * d
            v -= w
            m -= w * np.array([s["x_m"], s["y_m"], d / 2])
        return v, m / v


class BoxTank:
    def __init__(self, r, n=120):
        self.r = r
        self.name = r["tank"]
        self.rho = float(r["density_t_m3"])
        self.xa, self.xf = float(r["x_aft_m"]), float(r["x_fwd_m"])
        self.yp, self.ys = float(r["y_port_m"]), float(r["y_starboard_m"])
        self.zb, self.zt = float(r["z_bottom_m"]), float(r["z_top_m"])
        self.px, self.py = float(r["pipe_x_m"]), float(r["pipe_y_m"])
        self.s = float(r["sounding_m"])
        gx = np.linspace(self.xa, self.xf, n + 1)
        gy = np.linspace(self.yp, self.ys, n + 1)
        self.X, self.Y = np.meshgrid(0.5 * (gx[1:] + gx[:-1]), 0.5 * (gy[1:] + gy[:-1]), indexing="ij")
        self.dA = (self.xf - self.xa) * (self.ys - self.yp) / n / n

    def fill(self, a, slope_x, slope_y):
        return np.clip(a + slope_x * self.X + slope_y * self.Y - self.zb, 0.0, self.zt - self.zb)


def shared_level(boxes, vol, slope_x, slope_y):
    lo = min(t.zb for t in boxes) - 40
    hi = max(t.zt for t in boxes) + 40
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        v = sum(t.fill(mid, slope_x, slope_y).sum() * t.dA for t in boxes)
        lo, hi = (mid, hi) if v < vol else (lo, mid)
    a = 0.5 * (lo + hi)
    m, v = np.zeros(3), 0.0
    for t in boxes:
        h = t.fill(a, slope_x, slope_y)
        v += h.sum() * t.dA
        m += np.array([(h * t.X).sum(), (h * t.Y).sum(), (h * (t.zb + h / 2)).sum()]) * t.dA
    return v, m / v


def solve(folder, **changes):
    o = dict(SWITCHES, **changes)
    with open(os.path.join(folder, "ship.json")) as f:
        part = json.load(f)
    body = HullFromOffsets(folder, part, o)
    samples = np.array([float(r["density_kg_m3"]) for r in read_table(os.path.join(folder, "water_samples.csv"))])
    if o["sample_screen"] and samples.size >= 3:      # a sample far from the others is a spoilt bottle
        samples = samples[np.abs(samples - np.median(samples)) < 5.0]
    sea = samples.mean() / 1000.0
    lbp = float(part["lbp_m"])

    # waterplane through the mean mark readings (both sides, start and end of the day)
    scales = {m["mark"]: float(m["x_m"]) for m in part["draft_marks"]}
    if not o["mark_positions"]:
        xs_sorted = sorted(scales, key=scales.get)
        scales[xs_sorted[0]], scales[xs_sorted[-1]] = 0.0, lbp
        if len(xs_sorted) == 3:
            scales[xs_sorted[1]] = lbp / 2
    reads = {}
    for r in read_table(os.path.join(folder, "draft_readings.csv")):
        reads.setdefault((r["mark"], r["side"]), []).append((float(r["clock_s"]), float(r["reading_m"])))
    t_drafts = np.mean([t for v in reads.values() for t, _ in v])
    pts3 = []
    for (mk, side), v in reads.items():
        z = float(np.mean([r for _, r in v]))
        x = scales[mk]
        y = float(body.half_breadth(x, z)) * (1 if side == "starboard" else -1)
        pts3.append((x, y, z))
    pts3 = np.array(pts3)
    if o["mark_screen"]:
        # the taped freeboards are coarser but independent: a mark scale that disagrees with them is not used
        fb = read_table(os.path.join(folder, "freeboard_tapes.csv"))
        q = np.array([(float(r["x_m"]), float(body.half_breadth(float(r["x_m"]), body.D)) * (1 if r["side"] == "starboard" else -1),
                       body.D - float(r["freeboard_m"])) for r in fb])
        pf = np.linalg.lstsq(np.c_[np.ones(len(q)), q[:, 0], q[:, 1]], q[:, 2], rcond=None)[0]
        dev = pts3[:, 2] - (pf[0] + pf[1] * pts3[:, 0] + pf[2] * pts3[:, 1])
        good = np.abs(dev - np.median(dev)) < 0.06
        if good.sum() >= 4:
            pts3 = pts3[good]
    A = np.c_[np.ones(len(pts3)), pts3[:, 0], pts3[:, 1]]
    wl_a, slope_x, slope_y0 = np.linalg.lstsq(A, pts3[:, 2], rcond=None)[0]

    app = [(float(a["volume_m3"]), np.array([a["x_m"], a["y_m"], a["z_m"]], float)) for a in part["appendages"]]
    v_app = sum(v for v, _ in app)
    m_app = sum(v * c for v, c in app) if app else np.zeros(3)

    def float_at(a, b, c):
        v, cen = body.immersed(a, b, c)
        return sea * (v + v_app), (v * cen + m_app) / (v + v_app)

    delta_t, b0 = float_at(wl_a, slope_x, slope_y0)
    v_day, _ = body.immersed(wl_a, slope_x, slope_y0)

    # liquid_moments: the surface passes through the sounding at the pipe, parallel to the waterplane
    boxes = [BoxTank(r) for r in read_table(os.path.join(folder, "tank_soundings.csv")) if float(r["sounding_m"]) > 0]
    by = {t.name: t for t in boxes}
    pools = []
    seen = set()
    for t in boxes:
        if t.name in seen:
            continue
        pool = [t]
        m = re.search(r"cross-over valve to (\S+) open", t.r.get("remark", ""))
        if m and o["joined_tanks"] and m.group(1) in by:
            pool.append(by[m.group(1)])
        seen.update(g.name for g in pool)
        vol = 0.0
        for g in pool:
            if "pressed full" in g.r.get("remark", ""):
                vol += (g.zt - g.zb) * (g.xf - g.xa) * (g.ys - g.yp)
            elif o["pipe_plane"]:
                a_t = g.zb + g.s - slope_x * g.px - slope_y0 * g.py
                vol += g.fill(a_t, slope_x, slope_y0).sum() * g.dA
            else:
                vol += g.s * (g.xf - g.xa) * (g.ys - g.yp)
        pools.append((pool, vol))

    def liquid_moments(b, c):
        out = []
        for pool, vol in pools:
            rho_l = pool[0].rho
            if o["exact_liquid"]:
                v, cen = shared_level(pool, vol, b, c)
            else:   # upright centroid, then the free-surface shift i*tan(heel)/V of each surface
                v, cen = shared_level(pool, vol, b, 0.0)
                if o["liquid_moves"] and len(pool) == 1:
                    g = pool[0]
                    full = (g.zt - g.zb) * (g.xf - g.xa) * (g.ys - g.yp)
                    if vol < full * 0.999:
                        i_t = (g.xf - g.xa) * (g.ys - g.yp) ** 3 / 12
                        cen = cen + np.array([0.0, i_t * c / v, 0.0])
                elif o["liquid_moves"]:
                    i_t = sum((g.xf - g.xa) * ((g.ys - g.yp) ** 3 / 12 + (g.ys - g.yp) * ((g.ys + g.yp) / 2 - cen[1]) ** 2)
                              for g in pool)
                    cen = cen + np.array([0.0, i_t * c / v, 0.0])
            out.append((rho_l * v, cen))
        return out

    # the rest of what was on board on the day
    taken_off, put_on, hung = [], [], []
    for r in read_table(os.path.join(folder, "lightship_list.csv")):
        m = float(r["mass_t"])
        at = np.array([r["x_m"], r["y_m"], r["z_m"]], float)
        if r["status"] == "fit":
            put_on.append((m, at))
        elif r["hung_from_z_m"].strip() and o["hang_point"]:
            hung.append((m, np.array([r["hung_from_x_m"], r["hung_from_y_m"], r["hung_from_z_m"]], float)))
        else:
            taken_off.append((m, at))
    wgt = {r["weight"]: r for r in read_table(os.path.join(folder, "test_weights.csv"))}
    cell = {r["weight"]: float(r["load_cell_t"]) for r in read_table(os.path.join(folder, "load_cell.csv"))}
    for k, w in wgt.items():     # the load cell weighed each weight as it came aboard; the paint is nominal
        w["mass_t"] = cell[k] if o["load_cell"] else float(w["painted_mass_t"])
    arrs = read_table(os.path.join(folder, "arrangements.csv"))

    # heel of each condition: pendulum differences from condition 0, anchored on the marks' list
    batten = read_table(os.path.join(folder, "battens.csv"))
    clock = np.array([float(r["clock_s"]) for r in batten])
    lens = {p["pendulum"]: float(p["length_m"]) * 1000.0 for p in part["pendulums"]}
    avg = []
    for s in arrs:
        sel = (clock >= float(s["from_clock_s"])) & (clock <= float(s["to_clock_s"]))
        avg.append([np.mean([float(batten[i][p]) for i in np.flatnonzero(sel)]) / lens[p] for p in lens])
    avg = np.array(avg)
    inc = read_table(os.path.join(folder, "inclinometer_log.csv"))
    ti = np.array([float(r["clock_s"]) for r in inc])
    ang = np.array([float(r["heel_deg"]) for r in inc])
    incl_tan = np.array([np.tan(np.radians(np.mean(ang[(ti >= float(s["from_clock_s"])) & (ti <= float(s["to_clock_s"]))])))
                        for s in arrs])
    step_incl = incl_tan - incl_tan[0]
    step_pend = avg - avg[0]
    if o["inclinometer_only"]:
        delta = step_incl
    else:
        use = np.ones(step_pend.shape[1], bool)
        if o["pendulum_screen"]:
            # each pendulum must agree with the inclinometer (2 % class) on the size of the heels
            ratio = np.array([np.dot(step_pend[:, j], step_incl) / np.dot(step_incl, step_incl) for j in range(step_pend.shape[1])])
            use = np.abs(ratio - 1) < 0.05
            if not use.any():
                use = np.abs(ratio - np.median(ratio)) == np.abs(ratio - np.median(ratio)).min()
        delta = step_pend[:, use].mean(axis=1)
    centres = np.array([0.5 * (float(s["from_clock_s"]) + float(s["to_clock_s"])) for s in arrs])
    drift_rate = 0.0
    if o["drift"]:
        # the starting arrangement is read again later in the day: any change is drift, taken as linear in time
        same = [k for k, s in enumerate(arrs) if all(s[w] == arrs[0][w] for w in "ABCD")]
        if len(same) >= 2:
            drift_rate = np.polyfit(centres[same], delta[same], 1)[0]
    heel_tan = slope_y0 - drift_rate * (t_drafts - centres[0]) + delta - drift_rate * (centres - centres[0])

    lhs, rvec = [], []
    others_mass = sum(m for m, _ in taken_off) + sum(m for m, _ in hung) + sum(float(w["mass_t"]) for w in wgt.values())
    for k, s in enumerate(arrs):
        c = heel_tan[k]
        if o["exact_heel"]:
            a = wl_a
            for _ in range(6):     # same displaced volume at this heel
                v, _ = body.immersed(a, slope_x, c)
                v2, _ = body.immersed(a + 0.002, slope_x, c)
                a -= (v - v_day) / ((v2 - v) / 0.002)
            _, bk = float_at(a, slope_x, c)
        else:
            km = upright_km(body, wl_a, slope_x, slope_y0, v_day, float_at)
            bk = b0 + np.array([0.0, (km - b0[2]) * (c - slope_y0), 0.0])
        mom = np.zeros(3)
        mass = 0.0
        for m, at in taken_off + hung:
            mom += m * at
            mass += m
        for w in wgt.values():
            y = float(w["y_port_m"]) if s[w["weight"]] == "P" else float(w["y_starboard_m"])
            mom += float(w["mass_t"]) * np.array([float(w["x_m"]), y, float(w["z_m"])])
            mass += float(w["mass_t"])
        for m, cen in liquid_moments(slope_x, c if o["liquid_moves"] else slope_y0):
            mom += m * cen
            mass += m
        m_rest = delta_t - mass
        lhs.append([m_rest, m_rest * c])
        rvec.append(delta_t * (bk[1] + c * bk[2]) - (mom[1] + c * mom[2]))
        if k == 0:
            mom_0, m_rest0 = mom, m_rest
    y_rest, z_rest = np.linalg.lstsq(np.array(lhs), np.array(rvec), rcond=None)[0]
    arm = slope_x if o["trim_lever"] else 0.0
    x_rest = (delta_t * (b0[0] + arm * b0[2]) - (mom_0[0] + arm * mom_0[2])) / m_rest0 - arm * z_rest

    m_light = m_rest0 + sum(m for m, _ in put_on)
    mom_light = m_rest0 * np.array([x_rest, y_rest, z_rest]) + sum(m * at for m, at in put_on)
    g = mom_light / m_light
    return dict(light_mass_t=round(float(m_light), 2), light_lcg_m=round(float(g[0]), 4), light_vcg_m=round(float(g[2]), 4),
                tcg_m=round(float(g[1]), 4), test_displacement_t=round(float(delta_t), 2),
                heel_deg=[round(float(np.degrees(np.arctan(t))), 3) for t in heel_tan])


def upright_km(body, a, slope_x, slope_y, vol, float_at):
    """upright KM from a small heel (used only when exact heel is switched off)"""
    d = 0.002
    aa = a
    for _ in range(6):
        v, _ = body.immersed(aa, slope_x, slope_y + d)
        v2, _ = body.immersed(aa + 0.002, slope_x, slope_y + d)
        aa -= (v - vol) / ((v2 - v) / 0.002)
    _, b1 = float_at(aa, slope_x, slope_y + d)
    _, b0 = float_at(a, slope_x, slope_y)
    return b0[2] + (b1[1] - b0[1]) / d


if __name__ == "__main__":
    ans = solve(sys.argv[1])
    with open(sys.argv[2], "w") as f:
        json.dump(ans, f)
