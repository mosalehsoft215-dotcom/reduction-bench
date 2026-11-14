"""Murphree vapour efficiencies of the rectifying and stripping trays from one column test run.

    python3 tray_efficiency.py <run folder> <answer.json>

Reads the column sheet, the shift log (test points, their set points and windows, the remarks), the DCS trend
(window means) and the laboratory results (converted to mole fractions), rebuilds every test point with the
rigorous stage model in stagewise.py and fits the two efficiencies to all laboratory compositions, the bottoms flow
and the steam flow at once.

The flows named in a point entry are held by flow controllers on their meters, so a stream's true flow is its set
point divided by that meter's as-found factor (a post-run check in the remarks; 1 otherwise). A volumetric reflux
set point becomes a mass flow through the density of the reflux at its true temperature. With the reflux on drum
level and the steam named instead, each point's reflux is fitted so that the column takes up the heat of that
steam. A bare-tray loss that was not measured is fitted with the steam meter's help; a feed heated against the
bottoms gets its temperature from the exchanger's heat balance. Compositions follow every move as one first-order
lag; its time constant comes from the column thermocouples, or, when there are none, from the samples themselves.
"""
import csv
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
SETTLED_MIN = 70
from stagewise import Mixture, solve_column  # noqa: E402


def mixture_data(sheet):
    comp = sheet["system"]["components"]
    act = sheet["system"]["liquid_activity"]
    return {"light": comp["light"], "heavy": comp["heavy"],
            "nrtl": {k: act[k] for k in ("a12", "b12_K", "a21", "b21_K", "alpha")}}


def density(d, x1, T):
    v = 0.0
    for c, x in ((d["light"], x1), (d["heavy"], 1 - x1)):
        rho = c["density_20C_kg_m3"] * (1 - c["expansion_1_K"] * (T - 293.15))
        v += x * c["molar_mass_g_mol"] / 1000 / rho
    return (x1 * d["light"]["molar_mass_g_mol"] + (1 - x1) * d["heavy"]["molar_mass_g_mol"]) / 1000 / v


POINT = re.compile(r"Point (\d+): feed ([0-9.]+) kg/h on tray (\d+), (reflux|steam) ([0-9.]+) (kg/h|m3/h), "
                   r"distillate ([0-9.]+) kg/h(?:, side draw ([0-9.]+) kg/h)?\.")


def read_points(folder, sheet):
    log = [l for l in (folder / "shift_log.txt").read_text().splitlines() if l.strip()]
    pts = {}
    for line in log:
        stamp, text = line.split(" | ", 1)
        t = datetime.strptime(stamp, "%d.%m.%Y %H:%M")
        m = POINT.match(text)
        if m:
            sp = {"feed": float(m.group(2)), m.group(4): float(m.group(5)), "dist": float(m.group(7))}
            if m.group(8):
                sp["side"] = float(m.group(8))
            pts[int(m.group(1))] = {"feed_tray": int(m.group(3)), "changed": t, "sp": sp,
                                    "reflux_unit": m.group(6) if m.group(4) == "reflux" else None}
        m = re.match(r"Point (\d+) lined out", text)
        if m:
            pts[int(m.group(1))]["start"] = t
        m = re.match(r"Point (\d+) test window closed", text)
        if m:
            pts[int(m.group(1))]["end"] = t
    with open(folder / "dcs_trend.csv") as fh:
        rows = list(csv.reader(fh))
    head, rows = rows[0], rows[1:]
    times = np.array([datetime.strptime(r[0], "%Y-%m-%d %H:%M") for r in rows])
    vals = np.array([[float(v) for v in r[1:]] for r in rows])
    col = {h: i for i, h in enumerate(head[1:])}
    with open(folder / "lab_results.csv") as fh:
        lab = list(csv.DictReader(fh))
    basis = "mass" if sheet["laboratory"]["result"].startswith("wt") else "molar"
    d = mixture_data(sheet)
    M1, M2 = d["light"]["molar_mass_g_mol"], d["heavy"]["molar_mass_g_mol"]
    key = [k for k in lab[0] if k.endswith("_pct")][0]
    for p in pts.values():
        # flows and pressure follow their controllers within minutes; average only where they have settled
        lo = max(p["start"], p["changed"] + timedelta(minutes=SETTLED_MIN))
        sel = (times >= lo) & (times <= p["end"])
        p["mean"] = {h: vals[sel, i].mean() for h, i in col.items()}
        p["lab"] = []
        for s in lab:
            ts = datetime.strptime(s["taken"], "%Y-%m-%d %H:%M")
            if p["start"] <= ts <= p["end"]:
                v = float(s[key]) / 100
                x = v / M1 / (v / M1 + (1 - v) / M2) if basis == "mass" else v
                p["lab"].append((s["point"], x, (ts - p["changed"]).total_seconds() / 60.0))
    return [pts[k] for k in sorted(pts)], head


def read_remarks(folder):
    """Shift-log entries that are not point entries, as (time, text)."""
    out = []
    for line in (folder / "shift_log.txt").read_text().splitlines():
        if " | " not in line:
            continue
        stamp, text = line.split(" | ", 1)
        if not re.match(r"Point \d+", text):
            out.append((datetime.strptime(stamp, "%d.%m.%Y %H:%M"), text))
    return out


def meter_factors(remarks, tags):
    """'as found it read X % above/below' in a remark on a meter: reading = factor * true flow."""
    fac = {}
    for _, text in remarks:
        m = re.search(r"\b(FI\d+)\b.*?as found it read ([0-9.]+) % (above|below)", text)
        if m:
            sign = 1.0 if m.group(3) == "above" else -1.0
            for key, tag in tags.items():
                if tag.split("_")[0] == m.group(1):
                    fac[key] = 1.0 + sign * float(m.group(2)) / 100.0
    return fac


COLLECTED = re.compile(r"\b(FI\d+)\b.*?settings of point (\d+): in ([0-9.]+) min at a steady indication of ([0-9.]+) "
                       r"(kg/h|m3/h) the \w+ took up ([0-9.]+) kg(?: of \w+)?(?: at ([0-9.]+) C)?")


def collection_factors(remarks, tags, points, d):
    """A collection test: the meter's steady indication against the mass a weighed tank or pot took up. A volumetric
    reflux meter is compared with that mass as volume at the tank temperature and the reflux composition (the
    distillate samples of that point: a total condenser)."""
    fac = {}
    for _, text in remarks:
        m = COLLECTED.search(text)
        if not m:
            continue
        kg_h = float(m.group(6)) / (float(m.group(3)) / 60.0)
        if m.group(5) == "m3/h":
            p = points[int(m.group(2)) - 1]
            x = float(np.mean([v for n, v, _ in p["lab"] if n == "distillate"]))
            true = kg_h / density(d, x, float(m.group(7)) + 273.15)
        else:
            true = kg_h
        for key, tag in tags.items():
            if tag.split("_")[0] == m.group(1):
                fac[key] = float(m.group(4)) / true
    return fac


def thermocouple_offsets(remarks):
    """'<tag> ... it read X K low|high': returns the correction to add to that tag's readings."""
    out = {}
    for _, text in remarks:
        m = re.search(r"\b(TI\d+_[a-z0-9]+_C)\b.*?read ([0-9.]+) K (low|high)", text)
        if m:
            out[m.group(1)] = float(m.group(2)) * (1 if m.group(3) == "low" else -1)
    return out


def lag_minutes(folder, points):
    """Time constant with which the column thermocouples follow each change of settings (first-order, one per run);
    None when the trend holds no column thermocouple."""
    with open(folder / "dcs_trend.csv") as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    times = np.array([datetime.strptime(r[0], "%Y-%m-%d %H:%M") for r in rows[1:]])
    tcs = [i for i, h in enumerate(head) if re.fullmatch(r"TI1\d\d_tray\d+_C|TI199_reboiler_C", h)]
    if not tcs:
        return None
    vals = np.array([[float(r[i]) for i in tcs] for r in rows[1:]])
    segs = []
    for k in range(1, len(points)):
        t0 = points[k]["changed"]
        t1 = points[k + 1]["changed"] if k + 1 < len(points) else times[-1] + timedelta(minutes=1)
        sel = (times >= t0) & (times < t1)
        m = np.array([(t - t0).total_seconds() / 60.0 for t in times[sel]])
        segs.append((m, vals[sel]))
    if not segs:
        return 14.0
    n_tc = len(tcs)

    def res(q):
        tau = np.exp(q[0])
        out = []
        for s, (m, v) in enumerate(segs):
            a = q[1 + 2 * n_tc * s: 1 + 2 * n_tc * s + n_tc]
            b = q[1 + 2 * n_tc * s + n_tc: 1 + 2 * n_tc * (s + 1)]
            out.append((b + (a - b) * np.exp(-m[:, None] / tau) - v).ravel())
        return np.concatenate(out) / 0.05
    q0 = [np.log(20.0)]
    for m, v in segs:
        q0 += list(v[0]) + list(v[-1])
    fit = least_squares(res, q0)
    return float(np.clip(np.exp(fit.x[0]), 2.0, 200.0))


def steam_heat(heating):
    """J one kg of heating steam gives up: latent heat, plus the sensible heat of a condensate leaving subcooled."""
    lam = float(re.search(r"latent heat ([0-9.]+) kJ/kg", heating).group(1)) * 1000
    m = re.search(r"condensate leaves at ([0-9.]+) C \(cp ([0-9.]+) kJ/kg K\)", heating)
    if m:
        Tsat = float(re.search(r"steam at ([0-9.]+) C", heating).group(1))
        lam += float(m.group(2)) * 1000 * (Tsat - float(m.group(1)))
    return lam


class Run:
    def __init__(self, folder):
        self.sheet = json.loads((folder / "column_sheet.json").read_text())
        self.d = mixture_data(self.sheet)
        self.mix = Mixture(self.d)
        s = self.sheet
        c = s["column"]
        N = c["trays"]
        self.bare = np.zeros(N)
        self.UA_known = True
        ua = 0.0
        if isinstance(c["insulation"], dict):
            for j in c["insulation"]["bare_trays"]:
                self.bare[j - 1] = 1.0
            ua = c["insulation"]["UA_W_K_per_bare_tray"]
            self.UA_known = ua is not None          # null: the bare trays lose an unknown, common UA
        self.UA = ua if self.UA_known else None
        self.col = dict(N=N, feed_tray=None, side_tray=c["side_draw"]["tray"] if c["side_draw"] else None,
                        condenser=s["condenser"]["type"], dP_tray_Pa=c["pressure_drop_per_tray_Pa"],
                        dP_condenser_Pa=s["condenser"]["pressure_drop_Pa"],
                        dP_reboiler_Pa=s["reboiler"]["pressure_drop_from_bottom_tray_Pa"])
        self.points, self.head = read_points(folder, s)
        m = s["meters"]
        self.tag = {k: v["tag"] for k, v in m.items() if v}     # a meter out of service is null
        self.ptag = s["pressure"]["tag"]
        self.punit = s["pressure"]["unit"].split()[0]
        self.p_drum = s["pressure"]["location"].startswith("reflux drum")
        self.ph = c.get("feed_preheater")
        self.q_steam = steam_heat(s["reboiler"]["heating"])
        self.kettle_loss = float(re.search(r"([0-9.]+) %", s["reboiler"]["kettle_heat_loss"]).group(1)) / 100
        self.on_steam = "steam" in self.points[0]["sp"]
        self.guess = {}
        ts = s["laboratory"].get("tray_samples", {})
        self.P_flash = ts["site_barometric_pressure_kPa"] * 1000.0 if "site_barometric_pressure_kPa" in ts else None
        # the feed comes from one tank: its composition is the mean of every feed sample of the run
        self.zF = float(np.mean([x for p in self.points for n, x, _ in p["lab"] if n == "feed"]))
        self.remarks = read_remarks(folder)
        self.fac = meter_factors(self.remarks, self.tag)
        self.fac.update(collection_factors(self.remarks, self.tag, self.points, self.d))
        for tag, corr in thermocouple_offsets(self.remarks).items():
            for p in self.points:
                if tag in p["mean"]:
                    p["mean"][tag] += corr
        self.tau = lag_minutes(folder, self.points)       # None: no column thermocouples, fitted to the samples
        self.moves = [(p["changed"] - self.points[0]["changed"]).total_seconds() / 60.0 for p in self.points]

    def true_flow(self, p, key):
        """A controlled stream: its meter is held at the set point, so true = set point / meter factor."""
        return p["sp"][key] / self.fac.get(key, 1.0)

    def operating(self, p):
        mn = p["mean"]
        P = mn[self.ptag] * (1000.0 if self.punit == "kPa" else 100.0)
        if self.p_drum:
            P += self.col["dP_condenser_Pa"]
        return dict(P_top=P, zF=self.zF, F_kg=self.true_flow(p, "feed"),
                    TF=mn["TI102_feed_C"] + 273.15 if "TI102_feed_C" in mn else None,
                    TR=mn["TI104_reflux_C"] + 273.15, D_kg=self.true_flow(p, "dist"),
                    S_kg=self.true_flow(p, "side") if "side" in p["sp"] else 0.0, T_amb=mn["TI001_ambient_C"] + 273.15)

    def solve(self, k, c, op, E):
        """Stage model at efficiencies E. Near-total recovery points do not converge from a cold start at high
        efficiency, so on failure walk the efficiencies up from 0.5 and carry the profile along."""
        def attempt(Ea, guess):
            try:
                r = solve_column(self.mix, c, op, Ea[0], Ea[1], guess=guess)
            except (ValueError, FloatingPointError):
                return None
            return r if r["ok"] else None

        for guess in [g for g in (self.guess.get(k), self.guess.get(k - 1)) if g is not None] + [None]:
            r = attempt(E, guess)
            if r is not None:
                return r
        start = attempt((0.5, 0.5), None)
        if start is None:
            return None
        s, step, z = 0.0, 0.25, start["z"]
        E0 = np.array([0.5, 0.5])
        while s < 1.0 and step > 1e-3:
            t = min(1.0, s + step)
            r = attempt(E0 + t * (np.asarray(E) - E0), z)
            if r is None:
                step /= 2
                continue
            s, z = t, r["z"]
            step = min(0.5, step * 1.5)
        return r if s >= 1.0 else None

    def model(self, k, p, E, kD=1.0, UA=None, R_kg=None):
        o = self.operating(p)
        o["D_kg"] *= kD
        mix = self.mix
        c = dict(self.col, feed_tray=p["feed_tray"], heat_loss_W_K=self.bare * (self.UA if UA is None else UA))
        F = o["F_kg"] / 3600 / mix.molar_mass(o["zF"])
        TF = o["TF"] if o["TF"] is not None else p.get("TF_last", 330.0)
        xD, x0, xs = 0.9, 0.9, 0.5
        r = None
        for _ in range(40):
            if R_kg is not None:
                L0_kg_s = R_kg / 3600
            elif p["reflux_unit"] == "m3/h":
                L0_kg_s = self.true_flow(p, "reflux") * density(self.d, x0, o["TR"]) / 3600
            else:
                L0_kg_s = self.true_flow(p, "reflux") / 3600
            op = dict(F=F, zF=o["zF"], hF=mix.hL(o["zF"], TF), L0=L0_kg_s / mix.molar_mass(x0), TR=o["TR"],
                      D=o["D_kg"] / 3600 / mix.molar_mass(xD), P_top=o["P_top"], T_amb=o["T_amb"],
                      S=o["S_kg"] / 3600 / mix.molar_mass(xs))
            r = self.solve(k, c, op, E)
            if r is None:
                raise RuntimeError(f"stage model did not converge at point {k + 1}")
            self.guess[k] = r["z"]
            new = [r["xD"], r["x0"], r["x"][c["side_tray"] - 1] if c["side_tray"] else 0.5]
            done = max(abs(new[0] - xD), abs(new[1] - x0), abs(new[2] - xs)) < 1e-11
            if self.ph:
                # feed/bottoms exchanger, no losses: F (h_F - h_in) = B (h_kettle - h_out)
                T_in = p["mean"][self.ph["feed_in"]] + 273.15
                T_out = p["mean"][self.ph["bottoms_out"]] + 273.15
                gain = r["B"] * (mix.hL(r["xB"], r["T"][-1]) - mix.hL(r["xB"], T_out)) / F
                TF_new = brentq(lambda T: mix.hL(o["zF"], T) - mix.hL(o["zF"], T_in) - gain, T_in - 50, T_in + 250)
                done = done and abs(TF_new - TF) < 1e-8
                TF = TF_new
            if done:
                break
            xD, x0, xs = new
        p["TF_last"] = TF
        return r

    def flashed(self, x, T):
        """Liquid left when tray liquid x at its bubble point T is let down to P_flash with no heat exchange."""
        mix, P = self.mix, self.P_flash
        x, T = float(x), float(T)
        if T <= float(mix.bubble(np.array(x), np.array(P))[0]):
            return x
        h_in = float(mix.hL(x, T))
        lo = mix.dew_liquid(x, P)[0] + 1e-12

        def gap(xl):
            Tl, yl = mix.bubble(np.array(xl), np.array(P))
            v = (x - xl) / (float(yl) - xl)
            return (1 - v) * float(mix.hL(xl, Tl)) + v * float(mix.hV(yl, Tl)) - h_in
        return brentq(gap, lo, x, xtol=1e-14, rtol=1e-13)

    def sample_state(self, r, name):
        """Liquid composition (and temperature, for tray valves) at a sample point in a steady state."""
        if name == "distillate":
            return r["xD"], None
        if name == "bottoms":
            return r["xB"], None
        if name == "side draw":
            return r["x"][self.col["side_tray"] - 1], None
        j = int(name.split()[1]) - 1
        return r["x"][j], r["T"][j]

    # parameter vector: E_rect, E_strip, kD, then UA (if not measured), tau (if no thermocouples), reflux of every
    # point (if the reflux is on drum level)
    def unpack(self, par):
        i = 3
        UA = None
        tau = self.tau
        if not self.UA_known:
            UA, i = par[i], i + 1
        if self.tau is None:
            tau, i = par[i], i + 1
        R = list(par[i:i + len(self.points)]) if self.on_steam else [None] * len(self.points)
        return par[:2], par[2], UA, tau, R

    def start(self):
        lo, hi, x0 = [0.2, 0.2, 0.97], [1.0, 1.0, 1.03], [0.65, 0.65, 1.0]
        if not self.UA_known:
            lo.append(0.0); hi.append(60.0); x0.append(4.0)
        if self.tau is None:
            lo.append(2.0); hi.append(200.0); x0.append(30.0)
        if self.on_steam:
            for p in self.points:
                D = self.true_flow(p, "dist")
                lo.append(0.2 * D); hi.append(6.0 * D); x0.append(1.5 * D)
        return x0, (lo, hi)

    def residuals(self, par):
        E, kD, UA, tau, R = self.unpack(par)
        out = [(kD - 1) / 0.002]
        states = [self.model(k, p, E, kD, UA, R[k]) for k, p in enumerate(self.points)]
        steam_meter = self.tag.get("steam")
        for k, (p, r) in enumerate(zip(self.points, states)):
            B_kg = r["B"] * self.mix.molar_mass(r["xB"]) * 3600
            B_meas = p["mean"][self.tag["bottoms"]] / self.fac.get("bottoms", 1.0)
            out.append((B_kg - B_meas) / (0.0015 * B_kg))
            steam_kg = r["Q_reb"] * (1 + self.kettle_loss) / self.q_steam * 3600
            if self.on_steam:
                out.append((steam_kg - self.true_flow(p, "steam")) / (0.001 * steam_kg))
            elif steam_meter and steam_meter in p["mean"]:
                out.append((steam_kg - p["mean"][steam_meter] / self.fac.get("steam", 1.0)) / (0.0015 * steam_kg))
            for name, x, age in p["lab"]:
                if name == "feed":
                    continue
                xm, Tm = self.sample_state(r, name)
                if k > 0:
                    # the state at this point's move, carried from the earlier points through the same lag
                    xs, Ts = self.sample_state(states[0], name)
                    for j in range(1, k):
                        w = np.exp(-(self.moves[j + 1] - self.moves[j]) / tau)
                        xj, Tj = self.sample_state(states[j], name)
                        xs = xj + (xs - xj) * w
                        Ts = None if Ts is None else Tj + (Ts - Tj) * w
                    w = np.exp(-age / tau)
                    xm = xm + (xs - xm) * w
                    Tm = None if Tm is None else Tm + (Ts - Tm) * w
                if Tm is not None and self.P_flash is not None:
                    xm = self.flashed(xm, Tm)
                sig = 0.0003 + 0.004 * min(x, 1 - x) + 0.001 * x       # the laboratory's stated repeatability
                out.append((xm - x) / sig)
        return np.array(out)


def main():
    folder, answer = Path(sys.argv[1]), Path(sys.argv[2])
    run = Run(folder)
    x0, bounds = run.start()
    fit = least_squares(run.residuals, x0, bounds=bounds, diff_step=1e-4, xtol=1e-10)
    E, kD, UA, tau, R = run.unpack(fit.x)
    out = {"emv_rectifying": float(E[0]), "emv_stripping": float(E[1]), "points": len(run.points),
           "chi2": float(np.sum(fit.fun ** 2)), "lag_minutes": float(tau)}
    if UA is not None:
        out["UA_W_K_per_bare_tray"] = float(UA)
    if run.on_steam:
        out["reflux_kg_h"] = [float(v) for v in R]
    answer.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
