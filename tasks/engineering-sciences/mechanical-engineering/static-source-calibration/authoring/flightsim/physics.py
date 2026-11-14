"""Air, airspeed and pressure relations shared by the campaign generator.

Units: SI throughout (Pa, K, m, m/s, kg) unless a name says otherwise.
"""
import numpy as np

R_AIR = 287.053          # J/(kg K), dry air
G0 = 9.80665             # m/s2
GAMMA = 1.4
CP = GAMMA * R_AIR / (GAMMA - 1.0)
P_SL, T_SL, LAPSE = 101325.0, 288.15, 0.0065
A_SL = np.sqrt(GAMMA * R_AIR * T_SL)
KNOT, FOOT = 0.514444, 0.3048


H_TROP = 11000.0
T_TROP = T_SL - LAPSE * H_TROP
P_TROP = P_SL * (T_TROP / T_SL) ** (G0 / (R_AIR * LAPSE))


def pressure_altitude(p):
    """ISA pressure altitude (m) of a static pressure (Pa): troposphere and the isothermal layer above 11 km"""
    p = np.asarray(p, float)
    low = T_SL / LAPSE * (1.0 - (p / P_SL) ** (R_AIR * LAPSE / G0))
    high = H_TROP - R_AIR * T_TROP / G0 * np.log(p / P_TROP)
    return np.where(p >= P_TROP, low, high)


def isa_pressure(h):
    h = np.asarray(h, float)
    low = P_SL * np.maximum(1.0 - LAPSE * h / T_SL, 1e-9) ** (G0 / (R_AIR * LAPSE))
    high = P_TROP * np.exp(-G0 * (h - H_TROP) / (R_AIR * T_TROP))
    return np.where(h <= H_TROP, low, high)


def cas_from_qc(qc):
    return A_SL * np.sqrt(5.0 * ((np.asarray(qc) / P_SL + 1.0) ** (2.0 / 7.0) - 1.0))


def qc_from_cas(v):
    return P_SL * ((1.0 + 0.2 * (np.asarray(v) / A_SL) ** 2) ** 3.5 - 1.0)


def mach_from(pt, ps):
    return np.sqrt(5.0 * ((np.asarray(pt) / np.asarray(ps)) ** (2.0 / 7.0) - 1.0))


def viscosity(t):
    """Sutherland's law for air (Pa s)"""
    return 1.458e-6 * t ** 1.5 / (t + 110.4)


def vapour_pressure(t):
    tc = t - 273.15
    return 610.94 * np.exp(17.625 * tc / (tc + 243.04))


class Atmosphere:
    """Real air of one test day: temperature profile (lapse rate with an optional low inversion), relative
    humidity, sea-level pressure; pressure by hypsometric integration with the virtual temperature."""

    def __init__(self, p_msl, t_msl, lapse=0.0065, inversion=None, rh=0.5, top=15500.0, dz=2.0, tropopause_k=None):
        self.z = np.arange(-100.0, top, dz)
        t = t_msl - lapse * self.z
        if tropopause_k is not None:        # isothermal above the tropopause
            t = np.maximum(t, tropopause_k)
        if inversion:                       # (base height m, depth m, temperature gain K)
            base, depth, gain = inversion
            t = t + gain * np.clip((self.z - base) / depth, 0.0, 1.0)
        self.t = t
        self.rh = rh
        e = rh * vapour_pressure(t)
        p = np.empty_like(self.z)
        i0 = int(np.argmin(np.abs(self.z)))
        p[i0] = p_msl
        tv = t / (1.0 - 0.378 * e / p_msl)
        for i in range(i0 + 1, len(self.z)):
            tv_i = 0.5 * (t[i - 1] + t[i]) / (1.0 - 0.378 * 0.5 * (e[i - 1] + e[i]) / p[i - 1])
            p[i] = p[i - 1] * np.exp(-G0 * dz / (R_AIR * tv_i))
        for i in range(i0 - 1, -1, -1):
            tv_i = 0.5 * (t[i + 1] + t[i]) / (1.0 - 0.378 * 0.5 * (e[i + 1] + e[i]) / p[i + 1])
            p[i] = p[i + 1] * np.exp(G0 * dz / (R_AIR * tv_i))
        self.p = p
        self.e = e

    def vapour(self, z):
        return np.interp(z, self.z, self.e)

    def at(self, z):
        return np.interp(z, self.z, self.p), np.interp(z, self.z, self.t)


def source_error(cl, mach, coef):
    """static source error as a fraction of the true impact pressure: (ps_indicated - ps) / qc"""
    a, b, c = coef
    return a + b * cl + c * mach * mach


def solve_flight(p, t, cas_ind, weight, wing, coef):
    """true Mach number at which the static system with its error shows the indicated airspeed cas_ind
    (vectorised bisection); returns mach, total pressure, indicated static pressure, lift coefficient"""
    p, t, cas_ind, weight = np.broadcast_arrays(*map(np.asarray, (p, t, cas_ind, weight)))
    lo, hi = np.full(p.shape, 0.03), np.full(p.shape, 0.95)
    for _ in range(64):
        m = 0.5 * (lo + hi)
        pt = p * (1.0 + 0.2 * m * m) ** 3.5
        cl = weight * G0 / (0.5 * GAMMA * p * m * m * wing)
        ps_i = p + source_error(cl, m, coef) * (pt - p)
        low = cas_from_qc(pt - ps_i) < cas_ind
        lo, hi = np.where(low, m, lo), np.where(low, hi, m)
    m = 0.5 * (lo + hi)
    pt = p * (1.0 + 0.2 * m * m) ** 3.5
    cl = weight * G0 / (0.5 * GAMMA * p * m * m * wing)
    return m, pt, p + source_error(cl, m, coef) * (pt - p), cl


def corrections(ias_kt, hp_ft, weight, wing, coef):
    """position corrections for an indicated condition: (dVpc knots, dHpc feet), with Vc = Vi + dVpc and
    Hc = Hi + dHpc; indicated quantities from the static source, true ones from the free-stream pressure"""
    ps_i = isa_pressure(hp_ft * FOOT)
    pt = ps_i + qc_from_cas(ias_kt * KNOT)
    dp = 0.0
    for _ in range(200):
        p = ps_i - dp
        m = mach_from(pt, p)
        cl = weight * G0 / (0.5 * GAMMA * p * m * m * wing)
        new = source_error(cl, m, coef) * (pt - p)
        if abs(new - dp) < 1e-9:
            break
        dp = new
    p = ps_i - dp
    return float((cas_from_qc(pt - p) / KNOT) - ias_kt), float((pressure_altitude(p) - hp_ft * FOOT) / FOOT)
