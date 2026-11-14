"""Physics behind the M-50 export metering skid records.

Every number written into a gas-day folder comes from the true process through these functions:

  gas states      CoolProp 6.6.0, HEOS mixture backend (pure-fluid reference equations with the GERG-2008
                  departure functions), single-phase gas imposed; the analysis image ships the same build.
  flow profile    turbulent power law u = Uc (1 - r/R)^(1/n), n = 1.66 log10(Re), with a linear top-bottom
                  tilt (1 + a y/R) that leaves the bulk velocity unchanged; each chord sees its line average.
  ultrasonics     straight chordal paths at angle theta to the axis, axial flow,
                  t = L / (c -/+ v cos(theta)) + transducer delay.
  chromatograph   detector area = mol % x response factor x injected-amount factor.
"""
import numpy as np
import CoolProp.CoolProp as CP
from CoolProp import AbstractState

R_GAS = 8.314462618
P_BASE = 101325.0
T_BASE = 288.15

# peaks reported by GC-501, in report order
PEAKS = ["C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6+", "N2", "CO2"]
C6_SPLIT = {"n-Hexane": 0.47466, "n-Heptane": 0.35340, "n-Octane": 0.17194}
FLUIDS = ["Methane", "Ethane", "Propane", "IsoButane", "n-Butane", "Isopentane", "n-Pentane",
          "n-Hexane", "n-Heptane", "n-Octane", "Nitrogen", "CarbonDioxide"]

# ideal-gas molar mass (g/mol) and molar superior calorific value at 15 C (kJ/mol), as configured
ISO6976 = {
    "Methane": (16.04246, 891.51), "Ethane": (30.06904, 1562.06), "Propane": (44.09562, 2220.99),
    "IsoButane": (58.1222, 2870.45), "n-Butane": (58.1222, 2879.63), "Isopentane": (72.14878, 3531.52),
    "n-Pentane": (72.14878, 3538.43), "n-Hexane": (86.17536, 4198.06), "n-Heptane": (100.20194, 4857.00),
    "n-Octane": (114.22852, 5515.80), "Nitrogen": (28.0134, 0.0), "CarbonDioxide": (44.0095, 0.0),
}

CHORD_Y = np.array([0.809017, 0.309017, -0.309017, -0.809017])     # path 1 is the top chord
WEIGHTS = 0.4 * np.sin(np.radians([36.0, 72.0, 108.0, 144.0])) ** 2


def expand(peaks):
    """10 reported peaks (mol %) -> 12 mole fractions summing to one (C6+ split into nC6/nC7/nC8)."""
    p = dict(zip(PEAKS, peaks))
    x = [p["C1"], p["C2"], p["C3"], p["iC4"], p["nC4"], p["iC5"], p["nC5"],
         p["C6+"] * C6_SPLIT["n-Hexane"], p["C6+"] * C6_SPLIT["n-Heptane"], p["C6+"] * C6_SPLIT["n-Octane"],
         p["N2"], p["CO2"]]
    x = np.array(x, float)
    return x / x.sum()


def molar_hs(x12):
    """Ideal-gas molar superior calorific value at 15 C, kJ/mol."""
    return float(sum(xi * ISO6976[f][1] for xi, f in zip(x12, FLUIDS)))


def molar_mass(x12):
    return float(sum(xi * ISO6976[f][0] for xi, f in zip(x12, FLUIDS))) / 1000.0


class GasStates:
    """Thin wrapper around one CoolProp mixture state (composition set per call)."""

    def __init__(self):
        self.st = AbstractState("HEOS", "&".join(FLUIDS))
        self.st.specify_phase(CP.iphase_gas)

    def _set(self, x12):
        self.st.set_mole_fractions(list(np.clip(x12, 0.0, 1.0) / np.sum(x12)))

    def at_pt(self, x12, p_pa, t_k):
        """(molar density mol/m3, speed of sound m/s) at p, T."""
        self._set(x12)
        self.st.update(CP.PT_INPUTS, float(p_pa), float(t_k))
        return self.st.rhomolar(), self.st.speed_sound()

    def at_dt(self, x12, rho_mol, t_k):
        """(pressure Pa, speed of sound m/s) of a trapped gas of fixed molar density."""
        self._set(x12)
        self.st.update(CP.DmolarT_INPUTS, float(rho_mol), float(t_k))
        return self.st.p(), self.st.speed_sound()

    def base_density(self, x12):
        return self.at_pt(x12, P_BASE, T_BASE)[0]


def _chord_ratio(n, y, pts=100001):
    h = np.sqrt(1.0 - y * y)
    xs = np.linspace(0.0, h, pts)
    r = np.minimum(np.sqrt(xs * xs + y * y), 1.0)
    line_mean = np.trapezoid((1.0 - r) ** (1.0 / n), xs) / h
    area_mean = 2.0 * n * n / ((n + 1.0) * (2.0 * n + 1.0))
    return line_mean / area_mean


_N_GRID = np.arange(6.0, 16.01, 0.05)
_K_GRID = np.array([[_chord_ratio(n, y) for y in CHORD_Y] for n in _N_GRID])


def chord_ratios(reynolds, tilt):
    """Line-average velocity of each chord over the bulk velocity."""
    n = 1.66 * np.log10(max(reynolds, 1.0e4))
    k = np.array([np.interp(n, _N_GRID, _K_GRID[:, j]) for j in range(4)])
    return k * (1.0 + tilt * CHORD_Y)


def integration_factor(reynolds):
    """Weighted chord sum over bulk velocity for the symmetric profile (tilt cancels)."""
    return float(np.dot(WEIGHTS, chord_ratios(reynolds, 0.0)))


def transit_times(path_len, axial_len, c, v_line):
    """Upstream and downstream flight times (s) of one chord in the gas, no delays."""
    cos_t = axial_len / path_len
    return path_len / (c - v_line * cos_t), path_len / (c + v_line * cos_t)
