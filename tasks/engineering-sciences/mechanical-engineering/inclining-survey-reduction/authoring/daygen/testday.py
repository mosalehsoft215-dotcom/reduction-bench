"""One inclining day: the ship as found, the weight shifts, and every record the surveyor takes home."""
import numpy as np

from hull import Hull
from afloat import Ship

SHIFTS = [  # where the four inclining weights stand at each condition
    dict(A="P", B="P", C="S", D="S"),
    dict(A="S", B="P", C="S", D="S"),
    dict(A="S", B="S", C="S", D="S"),
    dict(A="P", B="P", C="S", D="S"),
    dict(A="P", B="P", C="P", D="S"),
    dict(A="P", B="P", C="P", D="P"),
    dict(A="P", B="P", C="S", D="S"),
]


def station_grid(h):
    p = h.p
    st = list(np.arange(0.0, p["lbp"] + 1e-6, p["lbp"] / 40))
    step = p["lbp"] / 40
    x = -step
    while x > p["stern_tip"]:
        st.insert(0, x)
        x -= step
    st.insert(0, p["stern_tip"])
    x = p["lbp"] + step
    while x < p["stem_tip"]:
        st.append(x)
        x += step
    st.append(p["stem_tip"])
    return np.round(np.array(st), 3)


def waterline_grid(h):
    D = h.D
    z = [0.0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5]
    v = 2.0
    while v < D - 1e-6:
        z.append(v)
        v += 0.5
    z.append(D)
    return np.round(np.array(z), 3)


class Day:
    def __init__(self, c):
        self.c = c
        self.rng = np.random.default_rng(c["seed"])
        self.hull = Hull(c["form"])
        self.ship = Ship(self.hull, c["tanks"], [(a["volume_m3"], (a["x"], a["y"], a["z"])) for a in c["appendages"]])

    # the ship as she floats: lightship less what is missing, plus everything on board
    def loads(self, shift):
        c = self.c
        ls = c["lightship"]
        m = ls["mass"]
        mom = np.array([ls["lcg"], ls["tcg"], ls["vcg"]]) * m
        for it in c["items"]:
            if it["status"] == "add":           # part of the lightship, not on board on the day
                m -= it["mass"]
                mom -= it["mass"] * np.array(it["at"])
        fixed = [(m, *(mom / m))]
        hanging = []
        for it in c["items"]:
            if it["status"] == "deduct":
                if it.get("hangs_from") is None:
                    fixed.append((it["mass"], *it["at"]))
                else:
                    hanging.append((it["mass"], it["hangs_from"], it["hangs_from"][2] - it["at"][2]))
        for w in c["weights"]:
            y = w["y_port"] if shift[w["id"]] == "P" else w["y_stbd"]
            fixed.append((w["mass"], w["x"], y, w["z"]))
        return fixed, hanging

    def tune_trim(self):
        """place the lightship LCG and TCG so that the ship floats at the catalogue's trim and list"""
        c = self.c
        vols = {t["name"]: t["volume"] for t in c["tanks"] if t["volume"] > 0}
        lbp = c["form"]["lbp"]
        ls = c["lightship"]

        def attitude():
            fixed, hanging = self.loads(SHIFTS[0])
            q = self.ship.float_(c["rho"], fixed, vols, hanging)
            return -q[1] * lbp, float(np.degrees(np.arctan(q[2])))

        def secant(key, which, target):
            x0, x1 = ls[key], ls[key] + (0.5 if key == "lcg" else 0.05)
            ls[key] = x0
            f0 = attitude()[which] - target
            ls[key] = x1
            f1 = attitude()[which] - target
            for _ in range(8):
                if abs(f1) < 1e-4 or f1 == f0:
                    break
                x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
                ls[key] = x2
                x0, f0, x1, f1 = x1, f1, x2, attitude()[which] - target
            ls[key] = round(x1, 3)

        for _ in range(2):
            secant("lcg", 0, c["trim_by_stern"])
            secant("tcg", 1, c.get("list_deg", 0.25))

    # ---- the day --------------------------------------------------------------------------
    def wind(self, t):
        """outside heeling moment (t.m, to starboard): a beam wind freshening steadily through the day"""
        return self.c.get("wind_tm_per_h", 0.0) * (t - self.c["clock"]["drafts_start"]) / 3600.0

    def schedule(self):
        c = self.c
        seg, t = [], c["clock"]["test_start"]
        for k in range(len(SHIFTS)):
            if k > 0:
                seg.append((t, t + c["clock"]["move_s"], k - 1, k, "move"))
                t += c["clock"]["move_s"]
            seg.append((t, t + c["clock"]["steady_s"], k, k, "steady"))
            t += c["clock"]["steady_s"]
        return seg, t

    def run(self):
        c = self.c
        self.tune_trim()
        vols = {t["name"]: t["volume"] for t in c["tanks"] if t["volume"] > 0}
        self.seg, self.t_end = self.schedule()
        mids = {k: 0.5 * (t0 + t1) for t0, t1, k, _, kind in self.seg if kind == "steady"}
        planes = []
        guess = None
        for k, sh in enumerate(SHIFTS):
            fixed, hanging = self.loads(sh)
            q = self.ship.float_(c["rho"], fixed, vols, hanging, guess, m_heel=self.wind(mids[k]))
            planes.append(q)
            guess = q
        self.planes = np.array(planes)
        self.mids = mids
        fixed, hanging = self.loads(SHIFTS[0])
        self.plane_start = self.ship.float_(c["rho"], fixed, vols, hanging, planes[0],
                                            m_heel=self.wind(c["clock"]["drafts_start"]))
        self.plane_end = self.ship.float_(c["rho"], fixed, vols, hanging, planes[0],
                                          m_heel=self.wind(c["clock"]["drafts_end"]))
        # heel per unit outside moment, for the slow change of heel within each reading interval
        q = self.ship.float_(c["rho"], fixed, vols, hanging, planes[0], m_heel=self.wind(mids[0]) + 50.0)
        self.dc_dm = (q[2] - planes[0][2]) / 50.0
        return self.planes

    def heel_series(self):
        """tan(heel) of the ship at 2 Hz through the experiment, weight moves included (no gusts)"""
        dt = 0.5
        times = np.arange(self.seg[0][0], self.t_end, dt)
        cvals = np.empty_like(times)
        moving = np.zeros(times.size, bool)
        for t0, t1, k0, k1, kind in self.seg:
            sel = (times >= t0) & (times < t1)
            if kind == "steady":
                cvals[sel] = self.planes[k0][2] + self.dc_dm * (self.wind(times[sel]) - self.wind(self.mids[k0]))
            else:
                u = np.clip((times[sel] - t0 - 60) / (t1 - t0 - 120), 0, 1)
                c0 = self.planes[k0][2] + self.dc_dm * (self.wind(times[sel]) - self.wind(self.mids[k0]))
                c1 = self.planes[k1][2] + self.dc_dm * (self.wind(times[sel]) - self.wind(self.mids[k1]))
                cvals[sel] = c0 + (c1 - c0) * u
                moving[sel] = True
        gust = np.zeros(times.size)
        a = np.exp(-dt / self.c["noise"]["gust_tau_s"])
        e = self.rng.normal(0, self.c["noise"]["gust_rad"] * np.sqrt(1 - a * a), times.size)
        for i in range(1, times.size):
            gust[i] = a * gust[i - 1] + e[i]
        log = [(k, t0, t1) for t0, t1, k, _, kind in self.seg if kind == "steady"]
        return times, cvals + gust, moving, log

    # ---- the records ----------------------------------------------------------------------
    def mark_readings(self, plane, clock):
        c, rng, h = self.c, self.rng, self.hull
        a, bx, cy = plane
        bad = c.get("faults", {}).get("mark_offset", {})
        rows = []
        for mk in c["marks"]:
            for side, sgn in (("port", -1), ("starboard", 1)):
                z = a + bx * mk["x"]
                for _ in range(3):
                    y = sgn * float(h.half_breadth(mk["x"], z))
                    z = a + bx * mk["x"] + cy * y
                # a scale painted out of place reads high or low by the same amount at every reading
                z += bad.get(f"{mk['name']}/{side}", 0.0)
                for k in range(c["draft_reads"]):
                    r = z + rng.normal(0, c["noise"]["draft_m"])
                    rows.append((clock + 40 * k + rng.integers(0, 20), mk["name"], side, round(r * 200) / 200))
        return rows

    def drafts(self):
        """marks read from a boat before the first move and after the last one"""
        rows = self.mark_readings(self.plane_start, self.c["clock"]["drafts_start"])
        rows += self.mark_readings(self.plane_end, self.c["clock"]["drafts_end"])
        rows.sort()
        return rows

    def freeboards(self):
        """deck edge to water, taped at a few stations on both sides right after the first draft reading"""
        c, rng, h = self.c, self.rng, self.hull
        a, bx, cy = self.plane_start
        rows = []
        clock = c["clock"]["drafts_start"] + 420
        for x in c["freeboard_x"]:
            for side, sgn in (("port", -1), ("starboard", 1)):
                y = sgn * float(h.half_breadth(x, h.D))
                fb = h.D - (a + bx * x + cy * y)
                rows.append((clock, x, side, round(fb + rng.normal(0, c["noise"]["freeboard_m"]), 2)))
                clock += 45
        return rows

    def dockwater(self):
        c, rng = self.c, self.rng
        rows = []
        spoilt = c.get("faults", {}).get("sample", {})
        for i, d in enumerate(c["samples_depth"]):
            v = spoilt.get(i + 1, c["rho"] * 1000)      # a bottle rinsed in fresh water reads near 1000
            rows.append((i + 1, d, round(v + rng.normal(0, 0.25), 1)))
        return rows

    def instruments(self):
        """pendulum battens (mm) and the electronic inclinometer (degrees) over the whole experiment"""
        c, rng = self.c, self.rng
        times, tanh, moving, log = self.heel_series()
        true_len = c.get("faults", {}).get("pendulum_length", {})
        pend = []
        for pd in c["pendulums"]:
            # the deflection follows the wire's real length, whatever length was written down
            L = true_len.get(pd["name"], pd["length_m"]) * 1000.0
            amp = c["noise"]["roll_mm"] * (1 + 4 * moving)
            phase = rng.uniform(0, 2 * np.pi)
            swing = amp * np.sin(2 * np.pi * times / c["roll_period_s"] + phase) * (1 + 0.3 * np.sin(2 * np.pi * times / 517.0))
            v = pd["zero_mm"] + L * tanh + swing + rng.normal(0, c["noise"]["read_mm"], times.size)
            pend.append(np.round(v, 1))
        inc = c["inclinometer"]
        sel = (np.round(times * 2) % 2) == 0                  # 1 Hz
        ang = np.degrees(np.arctan(tanh[sel])) * inc["scale"] + inc["zero_deg"]
        ang = ang + 0.02 * np.sin(2 * np.pi * times[sel] / c["roll_period_s"] + 1.0) * (1 + 4 * moving[sel])
        ang = np.round((ang + rng.normal(0, 0.004, ang.size)) / 0.01) * 0.01
        return times, pend, (times[sel], np.round(ang, 2)), log
