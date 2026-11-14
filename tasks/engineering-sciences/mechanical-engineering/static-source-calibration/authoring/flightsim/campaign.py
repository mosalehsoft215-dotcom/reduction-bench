"""Flight-test campaign generator: one aircraft, one test day, a list of test cards flown in order.

    python3 campaign.py <output folder> [names...]

Each card is recorded as one or more blocks (the recorder runs from shortly before a card's stable window to its
end). Inside a block the aircraft state is computed at 20 Hz and written at 2 Hz.

What the generator puts into the records:
- the static ports read the free-stream pressure plus (a + b CL + c M^2) times the impact pressure, with one set of
  coefficients per flap setting (CL from the actual weight, M the true Mach number);
- the pitot probe reads the free-stream total pressure;
- each pressure line delays its transducer like a first-order lag (laminar flow in a long thin tube, so the time
  constant goes as the gas viscosity over the line pressure): the ground-test value scaled by (runway-level
  pressure / line pressure); the lines run through the cabin, which keeps one temperature from the ground test on;
- transducers turn pressure into volts through their stated calibration; the static one reads high by its stated
  zero shift times (cabin temperature - 20 C);
- the total-temperature probe reads T (1 + 0.2 r M^2), r its recovery factor;
- the true airspeed is the Mach number times the speed of sound of the moist air, sqrt(1.4 R T_v) with the
  virtual temperature T_v; the aircraft's hygrometer records the dew point of the outside air;
- GPS gives the antenna's position and velocity; the antenna sits on the fuselage centreline, above and aft of the
  static ports in body axes, so its height above the ports is dz cos(theta) - dx sin(theta) at pitch theta; in level
  flight the pitch attitude is the angle of attack CL / CL_alpha + alpha0(flap);
- the pilot flies the indicated airspeed of the card (with small scatter) at constant height;
- tower fly-bys: the photo grid is graduated for the vertical plane of the runway centreline, so a pass off that
  plane reads (height of the static-port mark above the lens) x D / d, d the actual horizontal distance from the
  camera, measured square to the runway;
- the fuel totalizer reads high by its stated percentage.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

import physics as ph

HZ, OUT_EVERY = 20, 10          # simulation rate, written every 10th sample (2 Hz)
NOISE = dict(pressure_pa=1.0, tat_k=0.04, dew_k=0.3, pitch_deg=0.15, vel_ms=0.03, pos_m=0.2, height_m=0.25, grid_m=0.10, baro_pa=1.5,
             tower_t=0.1, volts=1e-5, ias_kt=0.4, heading_deg=0.3)


def wind_at(c, z, t_s):
    levels = np.array([lv[0] for lv in c["wind"]["levels"]])
    spd = np.interp(z, levels, [lv[1] for lv in c["wind"]["levels"]])
    frm = np.radians(np.interp(z, levels, [lv[2] for lv in c["wind"]["levels"]]))
    de, dn = c["wind"].get("drift_ms_per_h", (0.0, 0.0))
    hours = t_s / 3600.0
    return -spd * np.sin(frm) + de * hours, -spd * np.cos(frm) + dn * hours


class Recorder:
    def __init__(self, c, rng):
        self.c, self.rng = c, rng
        self.rows, self.tower = [], []
        self.t = 0.0
        self.fuel = 0.0
        self.air = ph.Atmosphere(**c["air"])

    # ------------------------------------------------------------------ fuel and weight
    def burn(self, seconds, cas_kt, z):
        f = self.c["aircraft"]["flow"]
        rate = f[0] + f[1] * (cas_kt - 150.0) - f[2] * z / 1000.0
        self.fuel += rate * seconds
        self.t += seconds

    def weight(self):
        a = self.c["aircraft"]
        return a["zero_fuel_kg"] + a["fuel_at_start_kg"] - self.fuel

    # ------------------------------------------------------------------ one block
    def block(self, label_mask, cas_kt, z, heading_deg=None, flyby=None, flap=0, label=""):
        """label_mask, cas_kt, z: arrays at HZ over the block; heading_deg for legs, flyby dict for passes"""
        c, rng = self.c, self.rng
        n = len(cas_kt)
        dt = 1.0 / HZ
        t = self.t + np.arange(n) * dt
        p, temp = self.air.at(z)
        w = np.empty(n)
        fuel0 = self.fuel
        f = c["aircraft"]["flow"]
        rate = f[0] + f[1] * (cas_kt - 150.0) - f[2] * z / 1000.0
        fuel = fuel0 + np.concatenate([[0.0], np.cumsum(rate[:-1] * dt)])
        w = c["aircraft"]["zero_fuel_kg"] + c["aircraft"]["fuel_at_start_kg"] - fuel
        ias = cas_kt + self.ou(n, NOISE["ias_kt"], 8.0)
        m, pt, ps_i, cl = ph.solve_flight(p, temp, ias * ph.KNOT, w, c["aircraft"]["wing_area_m2"], c["coef"][flap])
        e_air = self.air.vapour(z)
        t_virtual = temp / (1.0 - 0.378 * e_air / p)                 # speed of sound of the moist air
        tas = m * np.sqrt(ph.GAMMA * ph.R_AIR * t_virtual)
        # pitch attitude in level flight is the angle of attack: alpha = CL / CL_alpha + alpha0(flap)
        lift = c["aircraft"]["lift"]
        pitch = np.radians(np.degrees(cl / lift["cl_alpha_per_rad"]) + lift["alpha0_deg"][flap])
        lg = np.log(np.maximum(e_air, 1.0) / 610.94)
        dew = 243.04 * lg / (17.625 - lg)
        we, wn = wind_at(c, z, t)
        we = we + self.ou(n, c["wind"]["gust_ms"], 15.0)
        wn = wn + self.ou(n, c["wind"]["gust_ms"], 15.0)
        if flyby is None:
            psi = np.full(n, np.radians(heading_deg))
            ge, gn = tas * np.sin(psi) + we, tas * np.cos(psi) + wn
            east = np.cumsum(ge) * dt
            north = np.cumsum(gn) * dt
            east, north = east + self.pos[0], north + self.pos[1]
        else:
            rw = np.radians(c["site"]["runway_heading_deg"])
            d_along = np.array([np.sin(rw), np.cos(rw)])
            d_side = np.array([np.cos(rw), -np.sin(rw)])          # to the right of the runway direction
            lat = flyby["lateral_m"] + flyby["drift_ms"] * (np.arange(n) * dt - flyby["pass_s"])
            vy = flyby["drift_ms"]
            wv = np.stack([we, wn], 1)
            # along-track ground speed so that the air speed equals the TAS
            bb = -2.0 * (wv @ d_along)
            cc = (vy * vy) - 2.0 * vy * (wv @ d_side) + (wv * wv).sum(1) - tas * tas
            gs = (-bb + np.sqrt(bb * bb - 4.0 * cc)) / 2.0
            along = np.cumsum(gs) * dt
            along = along - along[int(flyby["pass_s"] * HZ)] + c["site"]["camera"]["along_m"]
            pos = np.outer(along, d_along) + np.outer(lat, d_side)
            east, north = pos[:, 0], pos[:, 1]
            ge = gs * d_along[0] + vy * d_side[0]
            gn = gs * d_along[1] + vy * d_side[1]
            psi = np.arctan2(ge - we, gn - wn)
            ip = int(flyby["pass_s"] * HZ)
            cam = c["site"]["camera"]
            dist = cam["distance_m"] + lat[ip]          # camera sits on the left of the runway direction
            grid = (z[ip] - cam["height_m"]) * cam["distance_m"] / dist + rng.normal(0, NOISE["grid_m"])
            pb, tb = self.air.at(c["site"]["barometer_height_m"])
            _, tt = self.air.at(c["site"]["runway_elevation_m"] + 2.0)
            self.tower.append(dict(card=label, time_s=round(float(t[ip]), 1),
                                   grid_m="" if c["site"].get("camera_failed") else round(float(grid), 2),
                                   baro_hpa=round(float(pb + rng.normal(0, NOISE["baro_pa"])) / 100.0, 2),
                                   temp_c=round(float(tt - 273.15 + rng.normal(0, NOISE["tower_t"])), 1),
                                   rh_pct=int(round(100 * c["air"]["rh"]))))
        # pressure lines: first-order lag, time constant scaled from the ground test by the line pressure
        lines = c["lines"]
        p_ground = self.air.at(c["site"]["runway_elevation_m"] + 2.0)[0]
        ptl = self.lag(pt, lines["total_s"] * p_ground, dt)
        psl = self.lag(ps_i, lines["static_s"] * p_ground, dt)
        tat = temp * (1.0 + 0.2 * c["recovery"] * m * m)
        self.pos = (east[-1], north[-1])
        for i in range(0, n, OUT_EVERY):
            self.rows.append(dict(
                time_s=round(float(t[i]), 2),
                east_m=round(float(east[i] + rng.normal(0, NOISE["pos_m"])), 1),
                north_m=round(float(north[i] + rng.normal(0, NOISE["pos_m"])), 1),
                gps_height_m=round(float(z[i] + c["gps_above_ports_m"] * np.cos(pitch[i])
                                         - c["gps_aft_of_ports_m"] * np.sin(pitch[i])
                                         + rng.normal(0, NOISE["height_m"])), 1),
                v_east=round(float(ge[i] + rng.normal(0, NOISE["vel_ms"])), 3),
                v_north=round(float(gn[i] + rng.normal(0, NOISE["vel_ms"])), 3),
                v_up=round(float(np.gradient(z, dt)[i] + rng.normal(0, NOISE["vel_ms"])), 3),
                heading_deg=round(float(np.degrees(psi[i]) % 360.0 + rng.normal(0, NOISE["heading_deg"])), 1),
                pitch_deg=round(float(np.degrees(pitch[i]) + rng.normal(0, NOISE["pitch_deg"])), 1),
                total_v=self.volts("total", ptl[i] + rng.normal(0, NOISE["pressure_pa"])),
                static_v=self.volts("static", psl[i] + rng.normal(0, NOISE["pressure_pa"])),
                tat_c=round(float(tat[i] - 273.15 + rng.normal(0, NOISE["tat_k"])), 2),
                dewpoint_c=round(float(dew[i] + rng.normal(0, NOISE["dew_k"])), 1),
                cabin_c=round(float(c["cabin_c"] + rng.normal(0, 0.05)), 1),
                fuel_used_kg=round(float(fuel[i] * (1.0 + c["aircraft"]["totalizer_high_pct"] / 100.0)), 1),
                flap_deg=flap,
                card=label if label_mask[i] else ""))
        self.fuel = fuel[-1] + rate[-1] * dt
        self.t = t[-1] + dt

    def ou(self, n, sigma, tau):
        x = np.empty(n)
        a = np.exp(-1.0 / (HZ * tau))
        x[0] = self.rng.normal(0, sigma)
        e = self.rng.normal(0, sigma * np.sqrt(1 - a * a), n)
        for i in range(1, n):
            x[i] = a * x[i - 1] + e[i]
        return x

    @staticmethod
    def lag(x, tau_times_p, dt):
        y = np.empty_like(x)
        y[0] = x[0]
        for i in range(1, len(x)):
            tau = tau_times_p / y[i - 1]
            y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1.0 - np.exp(-dt / tau))
        return y

    def volts(self, which, p):
        cal = self.c["transducers"][which]
        a0, a1, a2 = cal["kpa"]
        kpa = p / 1000.0
        if which == "static":
            kpa = kpa + cal["zero_shift_kpa_per_k"] * (self.c["cabin_c"] - 20.0)
        v = (-a1 + np.sqrt(a1 * a1 - 4 * a2 * (a0 - kpa))) / (2 * a2) if a2 else (kpa - a0) / a1
        return round(float(v + self.rng.normal(0, NOISE["volts"])), 5)


def fly(c, folder):
    rng = np.random.default_rng(c["seed"])
    rec = Recorder(c, rng)
    rec.pos = (0.0, 0.0)
    cards = []
    for k, (kind, spec) in enumerate(c["cards"], 1):
        rec.burn(spec.get("gap_min", 4.0) * 60.0, spec.get("ias", 180.0), 1500.0)
        flap = spec.get("flap", 0)
        if kind == "flyby":
            label = f"F{k:02d}"
            n = int(50 * HZ)
            z = np.full(n, c["site"]["runway_elevation_m"] + spec["height_m"]) + rec.ou(n, 0.25, 6.0)
            mask = np.zeros(n, bool)
            mask[int(26 * HZ):int(34 * HZ)] = True
            rec.block(mask, np.full(n, float(spec["ias"])), z, flyby=dict(lateral_m=spec["lateral_m"],
                      drift_ms=spec["drift_ms"], pass_s=30.0), flap=flap, label=label)
            cards.append(dict(card=label, method="tower fly-by", ias_kt=spec["ias"], pressure_altitude_ft="",
                              heading_deg=c["site"]["runway_heading_deg"], flap_deg=flap))
        elif kind == "fourleg":
            z0 = float(np.interp(ph.isa_pressure(spec["hp_ft"] * ph.FOOT), rec.air.p[::-1], rec.air.z[::-1]))
            for leg in range(4):
                label = f"L{k:02d}-{leg + 1}"
                hdg = (spec["heading_deg"] + 90.0 * leg + rng.normal(0, 2.0)) % 360.0
                n = int(110 * HZ)
                z = np.full(n, z0) + rec.ou(n, 1.0, 10.0)
                mask = np.zeros(n, bool)
                mask[int(30 * HZ):] = True
                rec.block(mask, np.full(n, float(spec["ias"])), z, heading_deg=hdg, flap=flap, label=label)
                rec.burn(50.0, spec["ias"], z0)
                cards.append(dict(card=label, method=f"four-leg GPS run, leg {leg + 1}", ias_kt=spec["ias"],
                                  pressure_altitude_ft=spec["hp_ft"], heading_deg=round(hdg), flap_deg=flap))
        elif kind == "accel":
            label = f"A{k:02d}"
            z0 = float(np.interp(ph.isa_pressure(spec["hp_ft"] * ph.FOOT), rec.air.p[::-1], rec.air.z[::-1]))
            run = (spec["v2"] - spec["v1"]) / spec["rate_kt_s"]
            n = int((25 + run + 5) * HZ)
            ts = np.arange(n) / HZ
            lo, hi = min(spec["v1"], spec["v2"]), max(spec["v1"], spec["v2"])
            cas = np.clip(spec["v1"] + spec["rate_kt_s"] * (ts - 25.0), lo, hi)
            z = np.full(n, z0) + rec.ou(n, 1.0, 10.0)
            mask = (ts >= 25.0) & (ts <= 25.0 + run)
            rec.block(mask, cas, z, heading_deg=spec["heading_deg"], flap=flap, label=label)
            kind_txt = "level acceleration" if spec["rate_kt_s"] > 0 else "level deceleration"
            cards.append(dict(card=label, method=kind_txt, ias_kt=f"{spec['v1']}-{spec['v2']}",
                              pressure_altitude_ft=spec["hp_ft"], heading_deg=spec["heading_deg"], flap_deg=flap))
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "flight.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rec.rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rec.rows)
    with open(folder / "cards.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["card", "method", "ias_kt", "pressure_altitude_ft", "heading_deg",
                                           "flap_deg"], lineterminator="\n")
        w.writeheader()
        w.writerows(cards)
    if rec.tower:
        with open(folder / "tower.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rec.tower[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rec.tower)
    truth = dict(recovery_factor=c["recovery"])
    for req in c["requested"]:
        dv, dh = ph.corrections(req["ias_kt"], req["pressure_altitude_ft"], req["weight_kg"],
                                c["aircraft"]["wing_area_m2"], c["coef"][req["flap_deg"]])
        truth[req["key"]] = dv if req["key"] == "dvpc_kt" else dh
    return truth
