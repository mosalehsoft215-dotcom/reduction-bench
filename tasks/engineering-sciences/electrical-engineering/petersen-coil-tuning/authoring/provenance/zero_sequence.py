"""Time-domain zero-sequence circuit of the compensated 20 kV network, solved exactly, for simulating the record.

All branches are to earth.
  - Feeder k: per-phase capacitances C_k,ph and the same leakage conductance on each phase. The phase voltages to
    earth are the stiff, zero-sum source voltages e_ph(t) plus the neutral displacement v0(t).
  - ASC1 sits between the neutral of the zig-zag earthing transformer and earth. From the network neutral its branch
    is R_s + L_s in series (R_s: coil winding plus R0/3 of the earthing transformer; L_s: L0/3 of it), then the coil
    inductance L_p. The resistor of the power auxiliary winding, while switched in, is across L_p (referred to the
    main winding). The measuring winding sees the voltage across L_p.
  - The fixed coils out on F4 and F9 are pure inductances from the star points of their own earthing transformers
    to earth, across the network capacitance whenever their feeders are closed (both are, in every state here).
  - The fault is a resistance from one phase to earth at the busbar.
Within each switching state the system is linear and time invariant; it is stepped with the matrix exponential of the
state matrix extended by one oscillator per source harmonic, so time is not discretised approximately.
"""
import numpy as np
from scipy.linalg import expm

A120 = np.exp(2j * np.pi / 3)
U_N = 20000.0 / np.sqrt(3.0)          # nominal phase-to-earth voltage, V
W50 = 2 * np.pi * 50.0


def source_phasors(p):
    """rms phasors of the three source voltages for every harmonic: {h: array(3)} (phase order L1, L2, L3)"""
    e1 = p["u_ll_day"] / np.sqrt(3.0) * np.exp(1j * np.radians(p["phase0_deg"]))
    neg = p["neg_seq"][0] * np.exp(1j * np.radians(p["neg_seq"][1]))
    pos_set = np.array([1.0, A120 ** 2, A120])
    neg_set = np.array([1.0, A120, A120 ** 2])
    out = {1: e1 * (pos_set + neg * neg_set)}
    for h, mag, deg, seq in p["harmonics"]:
        out[h] = e1 * mag * np.exp(1j * np.radians(deg)) * (pos_set if seq > 0 else neg_set)
    return out


def feeder_sums(p, names):
    """per-phase capacitance (F) and total leakage conductance (S) of a set of feeders"""
    c = np.zeros(3)
    g = 0.0
    for n in names:
        ic, iw, dev = p["feeders"][n]
        c += ic / (W50 * U_N) / 3.0 * (1.0 + np.array(dev))
        g += iw / U_N
    return c, g


class Circuit:
    """one switching state"""

    def __init__(self, p, names, fault, resistor):
        self.p, self.resistor = p, resistor
        self.c_ph, self.g = feeder_sums(p, names)
        self.c = self.c_ph.sum()
        self.gf = 1.0 / p["r_fault"] if fault else 0.0
        self.rs = p["coil_r"] + p["et_r0"] / 3.0
        self.ls = p["et_x0_50"] / W50 / 3.0
        self.lp, self.rr = p["coil_l"], p["resistor"]
        fixed = sum(a for f, a in p["fixed_coils"].items() if f in names)     # fixed coils leave with their feeder
        self.lf = U_N / (W50 * fixed) if fixed > 0 else np.inf
        self.harm = sorted(source_phasors(p))
        self.w = 2 * np.pi * p["freq"]
        self.matrix = self._matrix()

    def source_current(self, e):
        """rms phasor of the current the sources push into the neutral node at one harmonic"""
        return {h: np.sum((1j * h * self.w * self.c_ph + self.g / 3.0) * e[h]) + self.gf * e[h][self.p["fault_phase"]]
                for h in self.harm}

    def _matrix(self):
        """d/dt of [v0, i_n, i_L, i_fixed, cos, sin, ...]; without the resistor i_L follows i_n"""
        n = 4 + 2 * len(self.harm)
        a = np.zeros((n, n))
        j = self.source_current(source_phasors(self.p))
        a[0, 0] = -(self.g + self.gf) / self.c
        a[0, 1] = a[0, 3] = -1.0 / self.c
        a[3, 0] = 0.0 if np.isinf(self.lf) else 1.0 / self.lf
        for k, h in enumerate(self.harm):
            a[0, 4 + 2 * k] = -np.sqrt(2.0) * j[h].real / self.c
            a[0, 5 + 2 * k] = np.sqrt(2.0) * j[h].imag / self.c
            a[4 + 2 * k, 5 + 2 * k] = -h * self.w
            a[5 + 2 * k, 4 + 2 * k] = h * self.w
        if self.resistor:
            a[1, 0] = 1.0 / self.ls
            a[1, 1] = -(self.rs + self.rr) / self.ls
            a[1, 2] = self.rr / self.ls
            a[2, 1] = self.rr / self.lp
            a[2, 2] = -self.rr / self.lp
        else:
            a[1, 0] = 1.0 / (self.ls + self.lp)
            a[1, 1] = -self.rs / (self.ls + self.lp)
            a[2, 0], a[2, 1] = a[1, 0], a[1, 1]
        return a

    def measuring_winding(self, x):
        if self.resistor:
            return self.rr * (x[1] - x[2])
        return self.lp * (x[0] - self.rs * x[1]) / (self.ls + self.lp)


def simulate(p, t_grid):
    """states on the uniform grid t_grid (s); returns (measuring-winding voltage, coil-branch current) in primary units.
    Switching instants: fault at p['t_fault'], resistor at p['t_resistor'], feeder trip at p['t_trip']."""
    everything = list(p["feeders"])
    remaining = [n for n in everything if n != p["tripped_feeder"]]
    events = sorted([(p["t_fault"], "fault"), (p["t_resistor"], "resistor"), (p["t_trip"], "trip")])
    states, on = [(-np.inf, Circuit(p, everything, False, False))], set()
    for t_ev, what in events:
        on.add(what)
        names = remaining if "trip" in on else everything
        states.append((t_ev, Circuit(p, names, "fault" in on and "trip" not in on, "resistor" in on)))
    harm = states[0][1].harm
    w = states[0][1].w
    dt = t_grid[1] - t_grid[0]
    x = np.zeros(4 + 2 * len(harm))
    for k, h in enumerate(harm):
        x[4 + 2 * k], x[5 + 2 * k] = np.cos(h * w * t_grid[0]), np.sin(h * w * t_grid[0])
    um, cur = np.empty(len(t_grid)), np.empty(len(t_grid))
    idx = 0
    circuit = states[0][1]
    step = expm(circuit.matrix * dt)
    for i, t in enumerate(t_grid):
        um[i], cur[i] = circuit.measuring_winding(x), x[1]
        if i + 1 == len(t_grid):
            break
        t_next = t_grid[i + 1]
        while idx + 1 < len(states) and states[idx + 1][0] < t_next:
            t_ev, new = states[idx + 1]
            x = expm(circuit.matrix * (t_ev - t)) @ x
            if new.resistor and not circuit.resistor:
                x[2] = x[1]
            circuit, t, idx = new, t_ev, idx + 1
            step = None
        if step is None:
            x = expm(circuit.matrix * (t_next - t)) @ x
            step = expm(circuit.matrix * dt)
        else:
            x = step @ x
        if not circuit.resistor:
            x[2] = x[1]
    return um, cur


def line_voltages(p, t):
    """u_L1L2(t), u_L2L3(t) in primary volts"""
    e = source_phasors(p)
    w = 2 * np.pi * p["freq"]
    ab = sum(np.sqrt(2.0) * ((e[h][0] - e[h][1]) * np.exp(1j * h * w * t)).real for h in e)
    bc = sum(np.sqrt(2.0) * ((e[h][1] - e[h][2]) * np.exp(1j * h * w * t)).real for h in e)
    return ab, bc


def truth(p):
    """the graded quantities: normal switching state, resistor disconnected, 20 kV and 50 Hz"""
    c_ph, g = feeder_sums(p, list(p["feeders"]))
    c = c_ph.sum()
    rs = p["coil_r"] + p["et_r0"] / 3.0
    xs = p["et_x0_50"] / 3.0 + W50 * p["coil_l"]
    branch = 1.0 / (rs + 1j * xs)
    i_c = W50 * c * U_N
    i_l = -branch.imag * U_N
    i_w = (g + branch.real) * U_N
    i_fixed = sum(p["fixed_coils"].values())
    asym = 1j * W50 * np.sum(c_ph * U_N * np.array([1.0, A120 ** 2, A120]))
    y = g + 1j * W50 * c + branch - 1j * i_fixed / U_N
    return dict(capacitive_current_A=i_c, coil_current_A=i_l, detuning_percent=100.0 * (i_c - i_l - i_fixed) / i_c,
                active_current_A=i_w, neutral_displacement_V=abs(asym / y))
