"""Physics of a water-filled steel line under a hydrostatic test.

Everything the simulation uses lives here:
  - water as TEOS-10 seawater of the stated Absolute Salinity (gsw.rho_t_exact, in-situ temperature);
  - Lame thick-wall response of each pipe element to its local internal pressure and to its temperature,
    plane strain for buried pipe (held axially by the soil) and closed-end free pipe above ground;
  - hydrostatic pressure along the profile, integrated element by element with the local density;
  - free air at the line's summit, ideal and isothermal with the water around it.
"""
import numpy as np
import gsw

G = 9.80665              # m/s2
R_GAS = 8.314462618      # J/(mol K)
E_STEEL = 207.0e9        # Pa
NU_STEEL = 0.30
ALPHA_STEEL = 11.7e-6    # 1/K
T_REF = 20.0             # C, temperature at which the register's dimensions hold
P_STD = 101325.0         # Pa
T_STD = 288.15           # K


def density(sa, t_c, p_gauge_pa, p_baro_pa):
    """kg/m3 of the test water at in-situ temperature t_c and gauge pressure p_gauge_pa."""
    sea_dbar = (np.asarray(p_gauge_pa) + np.asarray(p_baro_pa) - P_STD) / 1.0e4
    return gsw.rho_t_exact(sa, t_c, sea_dbar)


def bore_volume(v0, a0, b0, buried, t_c, p_gauge_pa):
    """Internal volume (m3) of an element whose bore volume is v0 at T_REF and zero gauge pressure."""
    a2, b2 = a0 * a0, b0 * b0
    den = E_STEEL * (b2 - a2)
    dt = np.asarray(t_c) - T_REF
    p = np.asarray(p_gauge_pa)
    if buried:
        er = (1.0 + NU_STEEL) * p * ((1.0 - 2.0 * NU_STEEL) * a2 + b2) / den + ALPHA_STEEL * (1.0 + NU_STEEL) * dt
        ez = 0.0 * er
    else:
        er = p * ((1.0 - 2.0 * NU_STEEL) * a2 + (1.0 + NU_STEEL) * b2) / den + ALPHA_STEEL * dt
        ez = p * (1.0 - 2.0 * NU_STEEL) * a2 / den + ALPHA_STEEL * dt
    return v0 * (1.0 + er) ** 2 * (1.0 + ez)


class Line:
    """A test section cut into elements of at most `step` metres along the chainage."""

    def __init__(self, segments, profile, sa, step=50.0):
        ch, z = np.asarray(profile[0], float), np.asarray(profile[1], float)
        rows = []
        for k, s in enumerate(segments):
            n = max(1, int(np.ceil((s["to_m"] - s["from_m"]) / step)))
            edges = np.linspace(s["from_m"], s["to_m"], n + 1)
            for i in range(n):
                mid = 0.5 * (edges[i] + edges[i + 1])
                rows.append((k, edges[i], edges[i + 1], mid, float(np.interp(mid, ch, z))))
        self.seg = np.array([r[0] for r in rows])
        self.x0 = np.array([r[1] for r in rows])
        self.x1 = np.array([r[2] for r in rows])
        self.x = np.array([r[3] for r in rows])
        self.z = np.array([r[4] for r in rows])
        self.length = self.x1 - self.x0
        od = np.array([segments[k]["od_mm"] for k in self.seg]) / 1000.0
        wt = np.array([segments[k]["wt_mm"] for k in self.seg]) / 1000.0
        self.b0 = od / 2.0
        self.a0 = od / 2.0 - wt
        self.v0 = np.pi * self.a0 ** 2 * self.length
        self.buried = np.array([segments[k]["kind"] == "buried" for k in self.seg])
        self.sa = sa
        self.n = len(rows)

    def volumes(self, t_c, p_el):
        out = np.empty(self.n)
        for flag in (True, False):
            m = self.buried == flag
            if m.any():
                out[m] = bore_volume(self.v0[m], self.a0[m], self.b0[m], flag, t_c[m], p_el[m])
        return out

    def pressures(self, p_ref, z_ref, t_c, p_baro, iters=2):
        """Gauge pressure (Pa) at every element's centre when a gauge whose tapping is at elevation z_ref
        on the first element reads p_ref: the water column is walked along the chainage element by element
        with the local density."""
        p = p_ref + 1000.0 * G * (z_ref - self.z)
        for _ in range(iters):
            rho = density(self.sa, t_c, p, p_baro)
            dz = np.diff(self.z)
            rmid = 0.5 * (rho[1:] + rho[:-1])
            p0 = p_ref + rho[0] * G * (z_ref - self.z[0])
            p = p0 - np.concatenate([[0.0], np.cumsum(rmid * G * dz)])
        return p
