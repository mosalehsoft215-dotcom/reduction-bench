"""Forward model of a Thellier experiment on oriented specimens of one lava flow.

Frames. Geographic: x north, y east, z down. A core's specimen frame comes from its azimuth and plunge, and a tilted
flow first turns the whole rock about its strike. The oven's field always lies along the specimen's +z axis.

Remanence. A specimen's remanence sits in blocking-temperature bins with weights w (sum 1). A field B acting on a bin
leaves m0 w f(B) K b, where b is the field direction, K the specimen's TRM anisotropy tensor (mean eigenvalue 1) and
f(B) = Bs tanh(B / Bs) the field law; a very large Bs is the linear law.

Nature. The flow cooled in the ancient field over a time t_nat at the sample's depth, so the primary remanence carries
1 + c log10(t_nat / t_lab) with c the cooling-rate sensitivity. A viscous part lies along the present field and goes
below about 200 C. A later flow on top reheats the rock to a temperature, and the bins up to it are rewritten in the
field of that time.

Laboratory. A zero-field step to T erases every bin at or below T. An in-field step to T then writes those bins in the
oven's true field. A pTRM check repeats an in-field step at a lower T. Heating above t_alter raises the capacity of the
rock to take a lab remanence by (1 + rate (Tmax - t_alter) / 100); the primary remanence is untouched.
"""
import math

import numpy as np

ROOM = 20.0


def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_y(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def rot_x(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def unit_from(dec, inc):
    d, i = math.radians(dec), math.radians(inc)
    return np.array([math.cos(i) * math.cos(d), math.cos(i) * math.sin(d), math.sin(i)])


def core_matrix(azimuth, plunge):
    """Specimen coordinates to geographic: the core's z axis points down the drill hole, at this azimuth and plunge."""
    return rot_z(math.radians(azimuth)) @ rot_y(math.radians(90.0 - plunge))


def tilt_matrix(strike, dip):
    """Pre-tilt geographic to present geographic for a bed of this strike and dip (dip to the right of the strike)."""
    s = math.radians(strike)
    return rot_z(s) @ rot_x(math.radians(dip)) @ rot_z(-s)


def fabric(p_degree, plane_normal, rng, jitter=0.01):
    """Symmetric tensor, mean eigenvalue 1: the flattening plane's normal has the short axis, degree p_degree."""
    n = np.asarray(plane_normal, float)
    n = n / np.linalg.norm(n)
    short = 3.0 / (2.0 * p_degree + 1.0)
    long_ = p_degree * short
    k = long_ * np.eye(3) + (short - long_) * np.outer(n, n)
    noise = rng.normal(0.0, jitter, (3, 3))
    k = k + 0.5 * (noise + noise.T)
    return 3.0 * k / np.trace(k)


def spectrum(temps, centre_low, centre_high, share_low, spread, rng):
    """Bin weights over the step temperatures (last step takes what is left); two blocking populations."""
    edges = np.concatenate([[ROOM], temps])
    cdf = lambda t: (share_low / (1 + np.exp(-(t - centre_low) / spread))
                     + (1 - share_low) / (1 + np.exp(-(t - centre_high) / (0.45 * spread))))
    w = np.diff(cdf(edges))
    w[-1] += 1.0 - cdf(temps[-1])
    w = np.clip(w * (1.0 + rng.normal(0.0, 0.015, len(w))), 1e-5, None)
    return w / w.sum()


def field_law(b, bs):
    return bs * math.tanh(b / bs)


class Specimen:
    def __init__(self, name, m0, w, temps, k, to_geo, bs, primary, viscous, t_alter, rate, noise, floor, rng):
        self.name, self.m0, self.w, self.temps, self.k, self.to_geo, self.bs = name, m0, w, temps, k, to_geo, bs
        self.primary = primary                # vector per unit of w, specimen frame (A m2)
        self.viscous = viscous                # (vector, unblocking temperature)
        self.t_alter, self.rate = t_alter, rate
        self.noise, self.floor, self.rng = noise, floor, rng
        self.scale = float(np.linalg.norm(primary))
        self.lab = np.zeros(len(w))           # bins now carrying a lab remanence, as a capacity-weighted size
        self.t_max = ROOM
        self.oven_axis = np.array([0.0, 0.0, 1.0])
        self.rewritten = np.ones(len(w))      # size of each bin's natural remanence relative to the primary

    def _capacity(self):
        return 1.0 + self.rate * max(self.t_max - self.t_alter, 0.0) / 100.0

    def _read(self, removed_to):
        keep = self.temps > removed_to
        m = self.primary * float((self.w * self.rewritten)[keep].sum())
        if removed_to < self.viscous[1]:
            m = m + self.viscous[0]
        m = m + self.m0 * self.k @ self.oven_axis * float(self.lab.sum())
        return m + self.rng.normal(0.0, self.noise * self.scale + self.floor, 3)

    def zero_step(self, t):
        self.t_max = max(self.t_max, t)
        self.lab[self.temps <= t] = 0.0
        self.removed = max(getattr(self, "removed", ROOM - 1.0), t)
        return self._read(self.removed)

    def field_step(self, t, b_true):
        self.t_max = max(self.t_max, t)
        cover = self.temps <= t
        self.lab[:] = 0.0
        self.lab[cover] = self.w[cover] * field_law(b_true, self.bs) * self._capacity()
        self.removed = max(getattr(self, "removed", ROOM - 1.0), t)
        return self._read(self.removed)

    def full_trm(self, b_true, direction, factor=1.0):
        """a total TRM after the run, for the anisotropy, cooling and field tests (capacity as left by the run)."""
        m = self.m0 * self.k @ direction * field_law(b_true, self.bs) * self._capacity() * factor
        return m + self.rng.normal(0.0, self.noise * float(np.linalg.norm(m)) + self.floor, 3)
