"""An independent reduction of one station (written apart from the reference, sharing none of its code).

    python second_route.py <station folder> <answer.json>

- heights: the same spherical geometry; the range offset by grid search so that the 1-km binned median ln D of
  near-zenith echoes (< 22 deg) and of oblique ones (> 42 deg) lie on one curve;
- diffusion: medians of ln D in 1-km bins from 84 to 96 km (echoes under 25 dB), weighted quadratic through them (the slope from echoes within 35 deg of the zenith);
- temperature: absolute (pressure moved to 90 km geometric when given on the geopotential surface) and gradient;
  logged notes decide which route is in doubt, and the two are averaged otherwise;
- wind: in each 3 km slab and local-time hour, a two-component fit of the radial velocities; each mean wind is the
  plain average of its 24 hourly values, which removes the tides whatever the echo rate, and the 90 km tide
  amplitudes come from a 24 h + 12 h harmonic fit to the hourly zonal winds there.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

RADIUS = 6371.0
MG_OVER_R = 28.9e-3 * 9.5 / 8.314 * 1e3
DIFF_CONST = 6.39e-2 * 2.5e-4


def load(folder):
    f = Path(folder)
    sheet = json.loads((f / "station_sheet.json").read_text())
    raw = np.genfromtxt(f / "echoes.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    hours = np.array([int(s[11:13]) + int(s[14:16]) / 60 + int(s[17:19]) / 3600 for s in raw["time_utc"]])
    text = (f / "station_log.txt").read_text().lower() + " ".join(sheet["notes"]).lower()
    return sheet, raw, hours, text


def heights(rng_km, zen):
    c = np.cos(np.radians(zen))
    return np.sqrt(rng_km ** 2 + RADIUS ** 2 + 2 * rng_km * RADIUS * c) - RADIUS


def binned(z, y, lo=84, hi=96):
    edges = np.arange(lo, hi + 1)
    mid, med, cnt = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (z >= a) & (z < b)
        if s.sum() > 30:
            mid.append(0.5 * (a + b)); med.append(np.median(y[s])); cnt.append(s.sum())
    return np.array(mid), np.array(med), np.array(cnt)


def main(folder, out):
    sheet, E, lt, text = load(folder)
    lam = 299.792458 / sheet["radar"]["frequency_mhz"]
    lnD = np.log(lam ** 2 / (16 * math.pi ** 2 * E["decay_time_s"]))
    under = E["snr_db"] < 25.0
    near, far = E["zenith_deg"] < 22, E["zenith_deg"] > 42

    def mismatch(off):
        z = heights(E["range_km"] - off, E["zenith_deg"])
        m1, v1, _ = binned(z[under & near], lnD[under & near], 86, 94)
        m2, v2, _ = binned(z[under & far], lnD[under & far], 86, 94)
        common = np.intersect1d(m1, m2)
        return np.mean([(v1[m1 == c][0] - v2[m2 == c][0]) ** 2 for c in common])

    grid = np.arange(-4.0, 4.01, 0.05)
    off = grid[int(np.argmin([mismatch(g) for g in grid]))]
    off = off if abs(off) > 0.6 else 0.0
    z = heights(E["range_km"] - off, E["zenith_deg"])
    def profile(mask):
        mid, med, cnt = binned(z[mask], lnD[mask])
        x = mid - 90.0
        A = np.vstack([np.ones_like(x), x, x * x]).T * np.sqrt(cnt)[:, None]
        return np.linalg.lstsq(A, med * np.sqrt(cnt), rcond=None)[0]

    a0 = profile(under)[0]                                   # level at 90 km from every underdense echo
    a1 = profile(under & (E["zenith_deg"] < 35.0))[1]       # slope from near-zenith echoes, whose heights are sharpest
    ana = sheet["analysis_pressure"]
    p_given = np.mean(list(ana["daily_pa"].values()))
    t_abs = math.sqrt(math.exp(a0) * p_given / DIFF_CONST)
    if "geopotential" in ana["surface"]:
        dz = RADIUS * 90.0 / (RADIUS - 90.0) - 90.0
        for _ in range(6):
            t_abs = math.sqrt(math.exp(a0) * p_given * math.exp(MG_OVER_R * dz / t_abs) / DIFF_CONST)
    t_grad = (2 * sheet["climatological_temperature_gradient_90km_k_per_km"] + MG_OVER_R) / a1
    # wind: per 3 km slab, hourly two-component fits; the means of the 24 hourly values give the profiles, and a
    # 24 h + 12 h harmonic fit to the hourly zonal winds of the 90 km slab gives the tide amplitudes
    th, ph = np.radians(E["zenith_deg"]), np.radians(E["azimuth_deg"])
    if any(k in text for k in ("transmitter", "frequency", "configuration")) and abs(t_grad / t_abs - 1) > 0.015:
        temp, factor = t_grad, t_grad / t_abs        # the configured wavelength is in doubt
    elif "warming" in text:
        temp, factor = t_abs, 1.0                    # the climatological gradient is in doubt
    elif abs(t_grad / t_abs - 1) < 0.03:
        temp, factor = 0.5 * (t_abs + t_grad), 1.0
    else:
        temp, factor = t_abs, 1.0
    east, north, tide = {}, {}, {}
    for level in (84, 87, 90, 93, 96):
        band = np.abs(z - level) < 1.5
        hourly_u, hourly_v = [], []
        for h in range(24):
            s = band & (np.floor(lt) == h)
            X = np.vstack([np.sin(th[s]) * np.sin(ph[s]), np.sin(th[s]) * np.cos(ph[s]), z[s] - level]).T
            sol = np.linalg.lstsq(X, E["radial_velocity_m_s"][s], rcond=None)[0]
            hourly_u.append(sol[0])
            hourly_v.append(sol[1])
        east[str(level)] = float(np.mean(hourly_u)) * factor
        north[str(level)] = float(np.mean(hourly_v)) * factor
        if level == 90:
            hrs = np.arange(24) + 0.5
            A = np.column_stack([np.ones(24)] + [f(2 * math.pi * hrs / P) for P in (24.0, 12.0) for f in (np.cos, np.sin)])
            c = np.linalg.lstsq(A, np.array(hourly_u), rcond=None)[0]
            tide = {"diurnal": float(math.hypot(c[1], c[2])) * factor, "semidiurnal": float(math.hypot(c[3], c[4])) * factor}
    Path(out).write_text(json.dumps({"zonal_wind_m_s": east, "meridional_wind_m_s": north,
                                     "zonal_tide_amplitude_90km_m_s": tide, "temperature_90km_k": temp,
                                     "check": {"offset_km": float(off), "t_abs": t_abs, "t_grad": t_grad}}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
