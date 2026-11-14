"""Physical model of one visit with a dual-beam imaging polarimeter (half-wave plate + Wollaston prism).

Each source has sky-frame normalised Stokes (q, u). For an exposure with instrument position angle PA (east of
north of the detector +x axis) the instrument-frame Stokes are the sky ones rotated by -PA, plus the telescope's
instrumental polarisation: a constant on the optical axis and, off axis, a radial term k r^2 (r in units of
1000 px from the optical axis), both fixed in the instrument frame. The half-wave plate at angle a with zero-angle
offset phi0 modulates the two beams: o, e = I/2 (1 +/- eps [Q cos 4(a - phi0) + U sin 4(a - phi0)]) times the beam
transmissions, where eps is the modulation efficiency of the plate in that band. Aperture photometry at four radii
loses the light outside each aperture, by a Moffat profile whose width follows the seeing of the exposure; the
e beam's image is a little softer than the o beam's. Photon, sky and read noise are added.
"""
import numpy as np

FILTERS_UM = {"U": 0.365, "B": 0.44, "V": 0.55, "R": 0.66, "I": 0.80}
ZERO_DEG = {"U": 1.9, "B": -0.8, "V": 0.4, "R": 1.3, "I": 2.6}
APERTURES = (4, 7, 12, 20, 30)
BETA = 4.0
SCALE = 0.25          # arcsec per pixel
AXIS = (1024.0, 1024.0)


def enclosed(r, fwhm_px):
    a = fwhm_px / (2.0 * np.sqrt(2.0 ** (1.0 / BETA) - 1.0))
    return 1.0 - (1.0 + (r / a) ** 2) ** (1.0 - BETA)


SPLIT = -88.0         # the e image lies 88 px below (-y) the o image


def spill(d, r, fwhm_px):
    """Fraction of a Moffat image centred d px away that falls inside a circle of radius r."""
    if d > r + 6.0 * fwhm_px:
        return 0.0
    a = fwhm_px / (2.0 * np.sqrt(2.0 ** (1.0 / BETA) - 1.0))
    g = np.linspace(-r, r, 121)
    xx, yy = np.meshgrid(g, g)
    inside = xx ** 2 + yy ** 2 <= r * r
    prof = (BETA - 1.0) / (np.pi * a * a) * (1.0 + ((xx - d) ** 2 + yy ** 2) / (a * a)) ** (-BETA)
    return float(prof[inside].sum() * (g[1] - g[0]) ** 2)


def serkowski(p_max, lam_max, k, lam):
    return p_max * np.exp(-k * np.log(lam_max / lam) ** 2)


def to_instrument(q, u, pa_deg):
    c, s = np.cos(np.radians(2 * pa_deg)), np.sin(np.radians(2 * pa_deg))
    return q * c + u * s, -q * s + u * c


class Visit:
    def __init__(self, p, rng):
        self.p, self.rng = p, rng
        self.ip0 = rng.normal(0.0, 3e-4, 2)

    def sources(self):
        p, rng = self.p, self.rng
        pt, tt = p["target_p"], p["target_angle"]
        r, az = p["target_r"], p["target_az"]
        out = [dict(name="T", r=r, az=az, q=pt * np.cos(np.radians(2 * tt)), u=pt * np.sin(np.radians(2 * tt)), flux=4.0e7)]
        if p.get("neighbour"):                # a polarised star whose e image falls 14 px from the target's o image
            nb = p["neighbour"]
            out.append(dict(name="F00", r=r, az=az, dx=9.9, dy=9.9 - SPLIT, q=nb["p"] * np.cos(np.radians(2 * nb["angle"])),
                            u=nb["p"] * np.sin(np.radians(2 * nb["angle"])), flux=nb["share"] * 4.0e7))
        tx, ty = 1000.0 * r * np.cos(np.radians(az)), 1000.0 * r * np.sin(np.radians(az))
        for j in range(p["field_stars"]):
            while True:                      # the target field is chosen so that no other star's images are near the target's
                rr = np.sqrt(rng.uniform(0.02, 0.9))
                az_j = rng.uniform(0.0, 360.0)
                fx, fy = 1000.0 * rr * np.cos(np.radians(az_j)), 1000.0 * rr * np.sin(np.radians(az_j))
                if min(np.hypot(fx - tx, fy - ty - dy) for dy in (0.0, SPLIT, -SPLIT)) > 60.0:
                    break
            fac = rng.uniform(0.85, 1.15)
            th = p["isp_angle"] + rng.normal(0.0, 3.0)
            out.append(dict(name=f"F{j + 1:02d}", r=rr, az=az_j,
                            q=p["isp"] * fac * np.cos(np.radians(2 * th)) + rng.normal(0.0, 3e-4),
                            u=p["isp"] * fac * np.sin(np.radians(2 * th)) + rng.normal(0.0, 3e-4),
                            flux=10 ** rng.uniform(5.7, 6.8)))
        return out

    def observe(self, kind, sources, pa_set, stall, exptime, start, log, phot, first=0):
        """pa_set: rotator setting (logged). stall: degrees the field turned by the last exposure (0 normally)."""
        p, rng = self.p, self.rng
        angles = np.arange(16) * 22.5
        clock = start
        for i, a in enumerate(angles):
            exp = f"{kind[0].upper()}{first + i + 1:02d}"
            drift = stall * i / 15.0
            pa = pa_set + drift                                  # actual position angle of the detector +x axis
            seeing = max(0.55, p["seeing_mean"] + rng.normal(0.0, p["seeing_sd"]))
            sky = p["sky_per_px"] * exptime / 60.0
            transp = 1.0 - abs(rng.normal(0.0, 0.03))
            log.append(dict(exposure=exp, object=kind, hwp_deg=a, exptime_s=exptime, instrument_pa_deg=round(pa_set % 360.0, 2),
                            airmass=round(1.05 + 0.4 * rng.random(), 3), seeing_arcsec=round(seeing, 2),
                            sky_e_per_px=round(sky, 1), start_s=round(clock, 1)))
            clock += exptime + 35.0
            fw_o = seeing / SCALE
            fw_e = fw_o * p["e_softer"]
            beams = {}
            for s in sources:
                az = s["az"] - drift                              # the field turns on the detector
                x = AXIS[0] + 1000.0 * s["r"] * np.cos(np.radians(az)) + s.get("dx", 0.0)
                y = AXIS[1] + 1000.0 * s["r"] * np.sin(np.radians(az)) + s.get("dy", 0.0)
                rr = np.hypot(x - AXIS[0], y - AXIS[1]) / 1000.0
                ph = np.arctan2(y - AXIS[1], x - AXIS[0])
                qi, ui = to_instrument(s["q"], s["u"], pa)
                k = p["k_radial"] * rr ** 2
                qi += self.ip0[0] + k * np.cos(2 * ph)
                ui += self.ip0[1] + k * np.sin(2 * ph)
                m = p["efficiency"] * (qi * np.cos(np.radians(4 * (a - p["zero"]))) + ui * np.sin(np.radians(4 * (a - p["zero"]))))
                flux = s["flux"] * exptime / 60.0 * transp
                beams[s["name"]] = (x, y, 0.5 * flux * (1 + m) * 1.03, 0.5 * flux * (1 - m) * 0.97)
            for s in sources:
                x, y, fo, fe = beams[s["name"]]
                row = dict(exposure=exp, source=s["name"], x_px=round(x + rng.normal(0, 0.05), 2), y_px=round(y + rng.normal(0, 0.05), 2))
                for ap in APERTURES:
                    area = np.pi * ap * ap
                    o = fo * enclosed(ap, fw_o)
                    e = fe * enclosed(ap, fw_e)
                    for other, (x2, y2, fo2, fe2) in beams.items():     # light of other images inside these apertures
                        if other == s["name"]:
                            continue
                        o += fe2 * spill(np.hypot(x2 - x, y2 + SPLIT - y), ap, fw_e) + fo2 * spill(np.hypot(x2 - x, y2 - y), ap, fw_o)
                        e += fo2 * spill(np.hypot(x2 - x, y2 - SPLIT - y), ap, fw_o) + fe2 * spill(np.hypot(x2 - x, y2 - y), ap, fw_e)
                    row[f"o_r{ap}"] = round(o + rng.normal(0.0, np.sqrt(o + area * (sky + 25.0))), 1)
                    row[f"e_r{ap}"] = round(e + rng.normal(0.0, np.sqrt(e + area * (sky + 25.0))), 1)
                phot.append(row)
        return clock

    def run(self):
        p, rng = self.p, self.rng
        log, phot = [], []
        qx, ux = p["unpol_intrinsic"]          # the "unpolarised" standard's own sky polarisation (zero normally)
        unpol = [dict(name="S", r=0.0, az=0.0, q=qx, u=ux, flux=3.0e7)]
        pa_u = rng.uniform(0, 180)
        clock = self.observe("unpolarised_standard", unpol, pa_u, 0.0, 30.0, 600.0, log, phot)
        clock = self.observe("unpolarised_standard", unpol, pa_u + 45.0, 0.0, 30.0, clock + 300.0, log, phot, first=16)
        std = p["standard"]
        ps = serkowski(std["p_max"], std["lam_max"], std["k"], FILTERS_UM[p["filter"]])
        qs, us = ps * np.cos(np.radians(2 * std["angle"])), ps * np.sin(np.radians(2 * std["angle"]))
        clock = self.observe("polarised_standard", [dict(name="S", r=0.0, az=0.0, q=qs, u=us, flux=3.0e7)],
                             rng.uniform(0, 180), 0.0, 30.0, clock + 400.0, log, phot)
        field = self.sources()
        self.observe("target_field", field, p["pa_mid"], p["stall"], 60.0, clock + 600.0, log, phot)
        cat = []
        for s in field:                       # sky offsets from the field centre, from the acquisition catalogue
            ang = np.radians(s["az"] + p["pa_mid"])
            x0 = 1000.0 * s["r"] * np.cos(np.radians(s["az"])) + s.get("dx", 0.0)
            y0 = 1000.0 * s["r"] * np.sin(np.radians(s["az"])) + s.get("dy", 0.0)
            rho, phi = np.hypot(x0, y0) * SCALE, np.arctan2(y0, x0) + np.radians(p["pa_mid"])
            cat.append(dict(source=s["name"], east_arcsec=round(rho * np.sin(phi) + rng.normal(0, 0.02), 2),
                            north_arcsec=round(rho * np.cos(phi) + rng.normal(0, 0.02), 2)))
        return log, phot, cat
