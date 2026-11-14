"""Physics of the synthetic VSM records: shape factors, magnetisation curves, instrument and file writing.

Every sample is a uniaxial ferrimagnet (M-type hexaferrite) that is divided into +c and -c domains in zero field.
A field across c turns the magnetisation of all domains alike against the anisotropy energy
K1 sin^2(theta) + K2 sin^4(theta) and against the shape energy of the moment component along the field; the c
components of the domains cancel, so only the demagnetising factor along the field enters:

    H_applied = (H1 + N_f Ms) m + H2 m^3   for |m| < 1,   m = M_f / Ms,   H1 = 2 K1 / (mu0 Ms),   H2 = 4 K2 / (mu0 Ms)

and m = +-1 beyond H1 + H2 + N_f Ms. A field along c moves the domain walls, which are held by a pinning field
H_c, so the internal field stays at +-H_c until the sample is saturated: M = (H_applied -+ H_c) / N_c.

Demagnetising factors are the magnetometric (volume-averaged, uniform-magnetisation) ones: Aharoni's closed form for
a rectangular box, 1/3 for a sphere, and for a circular disc the axial factor from its two face charges,
N_axial = (2/t) * integral_0^inf J1(kR)^2 (1 - exp(-k t)) / k^2 dk, with the radial factor (1 - N_axial) / 2.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy import integrate, special

MU0 = 4e-7 * math.pi
OE = 1000.0 / (4.0 * math.pi)          # A/m per oersted


def box_factor(along, other1, other2):
    """Aharoni (1998): magnetometric factor of a rectangular box along the edge `along` (full edge lengths)."""
    c, a, b = along / 2.0, other1 / 2.0, other2 / 2.0
    r = math.sqrt(a * a + b * b + c * c)
    ab, bc, ac = math.hypot(a, b), math.hypot(b, c), math.hypot(a, c)
    s = ((b * b - c * c) / (2 * b * c)) * math.log((r - a) / (r + a))
    s += ((a * a - c * c) / (2 * a * c)) * math.log((r - b) / (r + b))
    s += (b / (2 * c)) * math.log((ab + a) / (ab - a)) + (a / (2 * c)) * math.log((ab + b) / (ab - b))
    s += (c / (2 * a)) * math.log((bc - b) / (bc + b)) + (c / (2 * b)) * math.log((ac - a) / (ac + a))
    s += 2 * math.atan(a * b / (c * r))
    s += (a ** 3 + b ** 3 - 2 * c ** 3) / (3 * a * b * c) + ((a * a + b * b - 2 * c * c) / (3 * a * b * c)) * r
    s += (c / (a * b)) * (ac + bc) - (ab ** 3 + bc ** 3 + ac ** 3) / (3 * a * b * c)
    return s / math.pi


def disc_axial_factor(diameter, thickness):
    """(2R/t) * integral of J1(x)^2 (1 - exp(-x t/R)) / x^2 over x = kR, Simpson on [0, 4000] plus the tail."""
    R, t = diameter / 2.0, thickness
    x = np.linspace(0.0, 4000.0, 800001)
    x[0] = 1e-12
    y = special.j1(x) ** 2 * (-np.expm1(-x * t / R)) / (x * x)
    body = integrate.simpson(y, x=x)
    return (2.0 * R / t) * (body + 1.0 / (2.0 * math.pi * 4000.0 ** 2))


def shape_factor(shape, axis):
    kind = shape["kind"]
    if kind == "sphere":
        return 1.0 / 3.0
    if kind == "box":
        e = shape["edges_mm"]
        others = [e[k] for k in "xyz" if k != axis]
        return box_factor(e[axis], *others)
    if kind == "disc":
        n_ax = disc_axial_factor(shape["diameter_mm"], shape["thickness_mm"])
        return n_ax if axis == shape["axis"] else 0.5 * (1.0 - n_ax)
    raise ValueError(kind)


def hard_axis_m(h, a, b):
    """Root of a m + b m^3 = h with |m| <= 1 (a > 0, b >= 0), vectorised."""
    h = np.asarray(h, float)
    sat = a + b
    m = np.clip(h / sat, -1.0, 1.0)
    inside = np.abs(h) < sat
    if b > 0:
        p, q = a / b, -h[inside] / b               # m^3 + p m + q = 0, one real root (p > 0)
        d = np.sqrt(q * q / 4.0 + p ** 3 / 27.0)
        m[inside] = np.cbrt(-q / 2.0 + d) + np.cbrt(-q / 2.0 - d)
    else:
        m[inside] = h[inside] / a
    return m


def easy_axis_m(h, n_c, ms, hc, falling):
    """Wall motion against a pinning field: M = (H -+ Hc)/N, clipped at +-Ms; returns M/Ms."""
    h = np.asarray(h, float)
    shift = np.where(falling, -hc, hc)
    return np.clip((h - shift) / (n_c * ms), -1.0, 1.0)


def sweep(h_max, fine_to, fine_step, coarse_step):
    """Field set points of one loop: +h_max -> -h_max -> +h_max, fine steps inside +-fine_to."""
    up = [0.0]
    while up[-1] < h_max - 1e-9:
        step = fine_step if up[-1] < fine_to else coarse_step
        up.append(min(h_max, up[-1] + step))
    half = np.array(up)
    down = np.concatenate([half[::-1], -half[1:]])           # +h_max .. -h_max
    back = -down[1:]                                         # -h_max .. +h_max
    h = np.concatenate([down, back])
    falling = np.concatenate([np.ones(down.size, bool), np.zeros(back.size, bool)])
    return h, falling


class Instrument:
    """Electromagnet VSM of one measuring day."""

    def __init__(self, rng, scale, field_noise=(80.0, 3e-4), moment_noise=(2e-9, 5e-4), holder_slope=-2.3e-14):
        self.rng = rng
        self.scale = scale                    # indicated moment / true moment on that day
        self.field_noise = field_noise        # A/m absolute, relative
        self.moment_noise = moment_noise      # A m^2 absolute, relative
        self.holder_slope = holder_slope      # A m^2 per A/m, the empty holder

    def read(self, h_true, m_true, temp_k):
        rng = self.rng
        sh = np.hypot(self.field_noise[0], self.field_noise[1] * np.abs(h_true))
        h_read = h_true + sh * rng.standard_normal(h_true.size)
        total = m_true + self.holder_slope * h_true
        sm = np.hypot(self.moment_noise[0], self.moment_noise[1] * np.abs(total))
        m_read = self.scale * total + sm * rng.standard_normal(h_true.size)
        t = temp_k + 0.06 * rng.standard_normal(h_true.size)
        return h_read, m_read, t


FIELD_UNITS = {"H_Oe": OE, "mu0H_mT": 1e-3 / MU0, "H_kA_per_m": 1e3}
MOMENT_UNITS = {"m_emu": 1e-3, "m_uAm2": 1e-6, "m_Am2": 1.0}


def write_csv(path, h, m, t, field_col, moment_col):
    hf, mf = FIELD_UNITS[field_col], MOMENT_UNITS[moment_col]
    lines = [f"{field_col},{moment_col},T_K"]
    for a, b, c in zip(h / hf, m / mf, t):
        lines.append(f"{a:.2f},{b:.7e},{c:.2f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=1) + "\n", encoding="utf-8", newline="\n")
