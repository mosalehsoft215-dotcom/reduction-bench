"""Reference reduction of one polarimetry visit: the target's degree and angle of linear polarisation.

    python3 pol_reduce.py <visit folder> <output file>

Steps
1. Photometry: each source is measured in the largest aperture that no other image reaches. The e beam's image is
   softer than the o beam's, so a small aperture loses a different share of the two beams and that share follows
   the seeing; but a large aperture can also take in the e image of a neighbour that sits one beam split away.
   The position angle of every target-field exposure comes from the stars on the detector against their catalogue
   offsets on the sky, not from the rotator setting.
2. Model: for every exposure i and source j, F = (o - e) / (o + e) = eps [Q cos 4(a - phi0) + U sin 4(a - phi0)] + c_j,
   where a is the plate angle, phi0 the band's zero-angle offset, c_j the source's beam-transmission term, and
   (Q, U) the instrument-frame Stokes: the sky Stokes rotated by -PA_i of that exposure, plus the instrumental
   polarisation. Everything below is linear least squares on that model, in "apparent" units (times eps).
3. Instrumental polarisation on axis: the unpolarised standard, observed at two rotator settings; the part fixed in
   the instrument frame is instrumental, the part that turns with the rotator belongs to the star.
4. Efficiency: the polarised standard's apparent degree over its catalogue degree in this band, the Serkowski law
   P(lam) = P_max exp(-K ln^2(lam_max / lam)) evaluated at the band's wavelength.
5. Off axis: all field stars at once. Each is modelled as one shared sky polarisation (interstellar, which is real
   and is kept in the target) plus an instrumental term k r^2 along the radius from the optical axis, fixed in the
   instrument frame. Only the instrumental term is removed from the target.
6. Target: its sky (q, u) with every instrumental term removed, divided by eps, gives P and the angle.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

AXIS_UNIT = 1000.0


def load(folder):
    sheet = json.loads((folder / "visit_sheet.json").read_text())
    with open(folder / "exposures.csv") as f:
        expo = {r["exposure"]: r for r in csv.DictReader(f)}
    with open(folder / "photometry.csv") as f:
        phot = list(csv.DictReader(f))
    return sheet, expo, phot


def serkowski(p_max, lam_max, k, lam):
    return p_max * np.exp(-k * np.log(lam_max / lam) ** 2)


class Visit:
    def __init__(self, folder):
        self.sheet, self.expo, phot = load(folder)
        with open(folder / "field_catalogue.csv") as f:
            self.cat = {r["source"]: (float(r["east_arcsec"]), float(r["north_arcsec"])) for r in csv.DictReader(f)}
        self.phi0 = self.sheet["hwp_zero_offset_deg"][self.sheet["filter"]]
        ax = self.sheet["optical_axis_px"]
        split = self.sheet["e_image_offset_px"]
        radii = sorted(self.sheet["aperture_radii_px"])
        worst = max(float(e["seeing_arcsec"]) for e in self.expo.values()) / self.sheet["pixel_scale_arcsec"]
        margin = 1.5 * 1.15 * worst
        # where each source's two images sit (mean over its exposures), per sequence
        where = {}
        for r in phot:
            key = (self.expo[r["exposure"]]["object"], r["source"])
            where.setdefault(key, []).append((float(r["x_px"]), float(r["y_px"])))
        where = {k: np.mean(v, axis=0) for k, v in where.items()}
        self.ap = {}
        for (obj, src), (x, y) in where.items():
            mine = [(x, y), (x + split[0], y + split[1])]
            others = [pt for (o2, s2), (x2, y2) in where.items() if o2 == obj and s2 != src
                      for pt in ((x2, y2), (x2 + split[0], y2 + split[1]))]
            gap = min([np.hypot(a[0] - b[0], a[1] - b[1]) for a in mine for b in others] + [1e9])
            clean = [rad for rad in radii if rad + margin <= gap]
            self.ap[(obj, src)] = clean[-1] if clean else radii[0]
        self.pa = self.measured_angles(phot, ax)
        self.phot, self.axis = phot, ax
        self.build_rows()

    def build_rows(self):
        ax = self.axis
        self.rows = []
        for r in self.phot:
            e = self.expo[r["exposure"]]
            ap = self.ap[(e["object"], r["source"])]
            o, eb = float(r[f"o_r{ap}"]), float(r[f"e_r{ap}"])
            dx, dy = (float(r["x_px"]) - ax[0]) / AXIS_UNIT, (float(r["y_px"]) - ax[1]) / AXIS_UNIT
            a = np.radians(4 * (float(e["hwp_deg"]) - self.phi0))
            pa = self.pa.get(r["exposure"], float(e["instrument_pa_deg"]))
            self.rows.append(dict(obj=e["object"], src=r["source"], F=(o - eb) / (o + eb), sig=2.0 / np.sqrt(o + eb),
                                  c4=np.cos(a), s4=np.sin(a), pa=np.radians(2 * pa),
                                  r2=dx * dx + dy * dy, ph=2 * np.arctan2(dy, dx)))

    def measured_angles(self, phot, ax):
        """Position angle of the detector +x axis in each target-field exposure, from the stars' positions against
        their catalogue offsets on the sky (sky angle = detector angle + PA). The rotator setting in the log says
        where it was told to be, the stars say where it was."""
        per = {}
        for r in phot:
            if self.expo[r["exposure"]]["object"] != "target_field" or r["source"] not in self.cat:
                continue
            dx, dy = float(r["x_px"]) - ax[0], float(r["y_px"]) - ax[1]
            if np.hypot(dx, dy) < 100.0:
                continue
            east, north = self.cat[r["source"]]
            d = np.arctan2(east, north) - np.arctan2(dy, dx)
            per.setdefault(r["exposure"], []).append(d)
        return {k: float(np.degrees(np.angle(np.mean(np.exp(1j * np.array(v)))))) for k, v in per.items()}

    def pick(self, obj, src=None, skip=None):
        return [r for r in self.rows if r["obj"] == obj and (src is None or r["src"] == src) and r["src"] != skip]

    @staticmethod
    def solve(rows, columns, known):
        """Weighted least squares of F - known(row) on columns(row); one transmission term per source."""
        names = sorted({r["src"] for r in rows})
        A, b, w = [], [], []
        for r in rows:
            off = [1.0 if r["src"] == n else 0.0 for n in names]
            A.append(columns(r) + off)
            b.append(r["F"] - known(r))
            w.append(1.0 / r["sig"])
        A, b, w = np.array(A), np.array(b), np.array(w)
        sol, *_ = np.linalg.lstsq(A * w[:, None], b * w, rcond=None)
        cov = np.linalg.pinv((A * w[:, None]).T @ (A * w[:, None]))
        return sol, cov

    @staticmethod
    def sky_columns(r):
        c, s = np.cos(r["pa"]), np.sin(r["pa"])
        # Q = q c + u s, U = -q s + u c  ->  F = Q c4 + U s4
        return [c * r["c4"] - s * r["s4"], s * r["c4"] + c * r["s4"]]

    def reduce(self):
        # 3. on-axis instrumental polarisation (apparent units), instrument frame
        # The standard was observed at two rotator settings. Instrumental polarisation stays fixed in the instrument
        # frame while any polarisation of the star itself turns with the rotator, so both are fitted together and
        # only the instrumental part is kept.
        un = self.pick("unpolarised_standard")
        turned = np.ptp([r["pa"] for r in un]) > 0.1
        cols = (lambda r: [r["c4"], r["s4"]] + self.sky_columns(r)) if turned else (lambda r: [r["c4"], r["s4"]])
        ip, _ = self.solve(un, cols, lambda r: 0.0)
        ip0 = ip[:2]
        base = lambda r: ip0[0] * r["c4"] + ip0[1] * r["s4"]
        # 4. efficiency from the polarised standard
        st, _ = self.solve(self.pick("polarised_standard"), self.sky_columns, base)
        cat = self.sheet["polarised_standard"]
        lam = self.sheet["filter_wavelength_um"][self.sheet["filter"]]
        eps = np.hypot(st[0], st[1]) / serkowski(cat["p_max_percent"] / 100.0, cat["lambda_max_um"], cat["k"], lam)
        # 5. field stars: shared sky polarisation + radial instrumental term
        radial = lambda r: [r["r2"] * (np.cos(r["ph"]) * r["c4"] + np.sin(r["ph"]) * r["s4"])]
        fs, _ = self.solve(self.pick("target_field", skip=self.sheet["target"]),
                           lambda r: self.sky_columns(r) + radial(r), base)
        k = fs[2]
        # 6. target
        known = lambda r: base(r) + k * radial(r)[0]
        tg, cov = self.solve(self.pick("target_field", self.sheet["target"]), self.sky_columns, known)
        q, u = tg[0] / eps, tg[1] / eps
        p2 = q * q + u * u - (cov[0, 0] + cov[1, 1]) / (2 * eps * eps)
        p = np.sqrt(max(p2, 0.0))
        return {"target_p_percent": float(100 * p), "target_angle_deg": float((0.5 * np.degrees(np.arctan2(u, q))) % 180.0),
                "efficiency": float(eps), "radial_term_percent": float(100 * k / eps),
                "target_aperture_px": self.ap[("target_field", self.sheet["target"])]}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(Visit(Path(sys.argv[1])).reduce(), indent=1) + "\n")
