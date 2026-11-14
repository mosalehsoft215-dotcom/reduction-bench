"""Reduces one vault calibration session of a moving-coil seismometer to its open-circuit generator constant and
the channel sensitivity at 5 Hz.

    python3 vault_reduce.py <session directory> <result.json>

1. Shape: both records are modelled whole, each as the superposition of exact step responses at the logged event
   times plus the digitizer's FIR delay, with a quadratic baseline per record.  The calibration-coil steps carry the
   damping of the coil terminated by the shunt and the digitizer input in parallel (h), the weight lifts that of the
   coil terminated by the digitizer input alone (h_o); the natural angular frequency w0 is shared.  One fit over w0,
   h and h_o, with the step amplitudes and the baselines solved linearly inside it.
2. Commands without a step: a first fit gives every logged event its own amplitude; an event whose amplitude is
   under half the median of its record's is a command the relay or the solenoid did not act on and is dropped (an
   event whose first 3 s lie mostly inside a disturbance left out by step 3 is kept).
3. Earthquake: residual samples above six times the noise level plus the median residual at the same delay after
   the other events (a shape misfit repeats at every event, a quake does not), widened by 4 s on either side, are
   left out and the fit is repeated (twice before and twice after step 2), so a quake inside a record does not bias
   it.
4. Amplitudes: lift scale L = (lift amplitude) / (m_w g) = S_d G Z / (R_c + Z) and coil-step scale
   C = (step amplitude) / I = S_d G R_p / (R_c + R_p) K, with Z the digitizer input resistance and
   R_p = R_s Z / (R_s + Z).
5. Both coils share one gap, so K = G l_cal / l_sig, which turns C / L into G as a function of R_c; the damping
   difference h - h_o = G^2 / (2 m w0) [1 / (R_c + R_p) - 1 / (R_c + Z)] then fixes R_c (its first root on a log grid
   from 10 ohm to 300 kohm, refined by Brent's method; with no root, the grid point closest to one) and with
   it G and S_d G = L (R_c + Z) / Z.  Neither the listed S_d, the sheet K nor the vault thermometer enters.
6. Sensitivity at 5 Hz = S_d G R_p / (R_c + R_p) |H(5 Hz)| with the fitted w0 and h.
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # small fits: one BLAS thread, no spin-waiting
import numpy as np  # noqa: E402
from scipy.optimize import brentq, least_squares  # noqa: E402

ALPHA_CU = 0.00393


def step_velocity(t, m, w0, h):
    wd = w0 * np.sqrt(1 - h * h)
    out = np.zeros_like(t)
    on = t > 0
    out[on] = np.exp(-h * w0 * t[on]) * np.sin(wd * t[on]) / (m * wd)
    return out


class Record:
    """One digitized record: its samples, the signed event onsets, the events kept and the samples used."""

    def __init__(self, counts, rate, onsets, signs):
        self.y = np.asarray(counts, float)
        self.t = np.arange(self.y.size) / rate
        self.onsets, self.signs = np.asarray(onsets), np.asarray(signs)
        self.keep = np.ones(self.onsets.size, bool)
        self.use = np.ones(self.y.size, bool)
        x = self.t / self.t[-1]
        self.base = np.column_stack([np.ones_like(x), x, x * x])

    def columns(self, m, w0, h):
        return np.column_stack([sg * step_velocity(self.t - t0, m, w0, h)
                                for t0, sg in zip(self.onsets, self.signs)])

    def solve(self, m, w0, h, each=False, full=False):
        """Amplitude (one shared by the kept events, or one per event) and the residual on the samples used (on
        every sample with full=True), the amplitudes and baseline fitted on the samples used."""
        cols = self.columns(m, w0, h)
        shape = cols if each else (cols[:, self.keep].sum(1, keepdims=True))
        a = np.column_stack([shape, self.base])
        coef = np.linalg.lstsq(a[self.use], self.y[self.use], rcond=None)[0]
        res = a @ coef - self.y
        return (coef[:shape.shape[1]] if each else float(coef[0])), (res if full else np.where(self.use, res, 0.0))


def fit(cal, lifts, m, start=None, each=False):
    """Joint fit: shared w0, terminated damping h (calibration steps) and h_o (lifts); returns w0, h, h_o.  With
    each=True every event has its own amplitude (a fit that a command without a step cannot bias)."""
    def res(p):
        w0, h, ho = p
        return np.r_[cal.solve(m, w0, h, each)[1], lifts.solve(m, w0, ho, each)[1]]
    bounds = ([0.3, 0.05, 0.001], [100.0, 0.98, 0.6])
    if start is not None:            # a refit from the previous solution
        return least_squares(res, start, bounds=bounds, x_scale=[1.0, 0.1, 0.01], xtol=1e-12, ftol=1e-12,
                             gtol=1e-12).x
    # starting w0 from the ringing of the lightly damped lifts: the median over the 16 s windows after each lift of
    # their spectra, so that a disturbance inside one window does not pull it
    dt = lifts.t[1] - lifts.t[0]
    n_win = int(round(14.0 / dt))
    specs = []
    for t0 in lifts.onsets:
        i0 = int(np.searchsorted(lifts.t, t0))
        seg = lifts.y[i0:i0 + n_win]
        if seg.size == n_win:
            specs.append(np.abs(np.fft.rfft((seg - seg.mean()) * np.hanning(n_win), 8 * n_win)))
    freq = np.fft.rfftfreq(8 * n_win, dt)
    spec = np.median(specs, axis=0)
    band = (freq > 0.4) & (freq < 3.0)
    w_start = 2 * np.pi * freq[band][np.argmax(spec[band])]
    best = None
    for h0 in (0.45, 0.75):
        for ho0 in (0.03, 0.12):
            sol = least_squares(res, [w_start, h0, ho0], bounds=bounds,
                                x_scale=[1.0, 0.1, 0.01], xtol=1e-12, ftol=1e-12, gtol=1e-12)
            if best is None or sol.cost < best.cost:
                best = sol
    return best.x


def drop_missed(rec, m, w0, h):
    """Drops the events whose own amplitude is under half the median of the record's; an event whose first 3 s
    are mostly left out (inside a disturbance) is kept, its amplitude being undetermined there."""
    amps = np.abs(rec.solve(m, w0, h, each=True)[0])
    seen = np.array([rec.use[(rec.t > t0) & (rec.t < t0 + 3.0)].mean() for t0 in rec.onsets])
    judged = seen > 0.5
    rec.keep = ~judged | (amps >= 0.5 * np.median(amps[judged]))


def mask_quake(rec, m, w0, h, each=False):
    """Leaves out a disturbance the event model cannot explain: a residual sample counts as disturbed when it exceeds
    six times the sum of the noise level and the median residual at the same delay after the other events (a misfit
    of the step shape repeats at every event, a quake does not); disturbed stretches are widened by 4 s."""
    res = np.abs(rec.solve(m, w0, h, each, full=True)[1])
    dt = rec.t[1] - rec.t[0]
    period = int(round(np.median(np.diff(rec.onsets)) / dt))
    first = int(np.ceil(rec.onsets[0] / dt))
    rows = [res[first + k * period:first + (k + 1) * period] for k in range(rec.onsets.size)]
    rows = np.array([np.r_[r, np.full(period - r.size, np.nan)] for r in rows])
    typical = np.nanmedian(rows, axis=0)
    noise = 1.4826 * np.median(res)
    limit = np.full(res.size, 6 * (noise + np.nanmedian(typical)))
    for k in range(rec.onsets.size):
        seg = limit[first + k * period:first + (k + 1) * period]
        seg[:] = 6 * (noise + typical[:seg.size])
    bad = res > limit
    width = int(round(4.0 / dt))
    if bad.any():
        bad = np.convolve(bad.astype(float), np.ones(2 * width + 1), mode="same") > 0
    rec.use = ~bad


def h5(w0, h, f=5.0):
    r = f / (w0 / (2 * np.pi))
    return r * r / np.sqrt((1 - r * r) ** 2 + (2 * h * r) ** 2)


def reduce_session(folder, missed=True, quake=True, zin=True, g_route="reciprocity", rc_route="joint", fir=True):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py."""
    folder = Path(folder)
    sensor = json.loads((folder / "sensor.json").read_text())
    sess = json.loads((folder / "session.json").read_text())
    with open(folder / "events.csv", newline="") as fh:
        ev = list(csv.DictReader(fh))
    rate = sess["sample_rate_hz"]
    delay = sess["digitizer_fir_delay_s"] if fir else 0.0
    m, rs = sensor["moving_mass_kg"], sensor["shunt_resistance_ohm"]
    z = sess["digitizer_input_impedance_ohm"] if zin else 1e15

    def record(name, fname):
        rows = [e for e in ev if e["record"] == name]
        onsets = [float(e["time_s"]) + delay for e in rows]
        signs = [1.0 if e["event"] in ("current on", "weight lifted") else -1.0 for e in rows]
        return Record(np.loadtxt(folder / fname, skiprows=1), rate, onsets, signs)

    cal, lifts = record("calibration", "calibration_record.csv"), record("open_circuit", "open_circuit_record.csv")
    # per-event amplitudes first: a command without a step then gets a zero amplitude instead of biasing the shape
    p = fit(cal, lifts, m, each=missed)
    for _ in range(2 if quake else 0):
        mask_quake(cal, m, p[0], p[1], missed)
        mask_quake(lifts, m, p[0], p[2], missed)
        p = fit(cal, lifts, m, p, each=missed)
    if missed:
        drop_missed(cal, m, p[0], p[1])
        drop_missed(lifts, m, p[0], p[2])
    p = fit(cal, lifts, m, p)
    for _ in range(2 if quake else 0):
        mask_quake(cal, m, p[0], p[1])
        mask_quake(lifts, m, p[0], p[2])
        p = fit(cal, lifts, m, p)
    w0, h, ho = p
    a_cal, a_lift = cal.solve(m, w0, h)[0], lifts.solve(m, w0, ho)[0]
    lift = a_lift / (sess["weight_lift_mass_g"] * 1e-3 * sess["local_gravity_m_per_s2"])   # S_d G Z / (R_c + Z)
    coil = a_cal / (sess["calibrator_current_ma"] * 1e-3)                                  # S_d G R_p / (R_c + R_p) K
    ratio = sensor["signal_coil_wire_in_gap_m"] / sensor["calibration_coil_wire_in_gap_m"]
    rp = rs * z / (rs + z)

    def g_of(rc):
        return ratio * coil / lift * z * (rc + rp) / (rp * (rc + z))

    def damping_gap(rc):
        return g_of(rc) ** 2 * (1 / (rc + rp) - 1 / (rc + z)) - 2 * m * w0 * (h - ho)

    rc20 = sensor["signal_coil_resistance_ohm_at_20c"]
    rc_vault = rc20 * (1 + ALPHA_CU * (sess["vault_temperature_c"] - 20.0))
    if rc_route == "joint":                               # first sign change on a log grid, then Brent's method
        grid = np.geomspace(10.0, 3e5, 120)
        vals = np.array([damping_gap(x) for x in grid])
        cross = np.flatnonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))
        if cross.size:
            rc = brentq(damping_gap, grid[cross[0]], grid[cross[0] + 1])
        else:                                             # no root: the grid point that comes closest
            rc = float(grid[np.argmin(np.abs(vals))])
    else:
        rc = rc_vault
    g = g_of(rc)
    if g_route == "lift":                                  # listed S_d taken as it stands
        g = lift * (rc + z) / z / sess["digitizer_counts_per_volt"]
    elif g_route == "sheet-k":                             # reciprocity with the sheet motor constant
        g = sensor["calibration_coil_motor_constant_n_per_a"] * ratio
    scale = lift * (rc + z) / z                            # S_d G
    return {"generator_constant_v_s_per_m": float(g),
            "sensitivity_5hz_counts_per_m_per_s": float(scale * rp / (rc + rp) * h5(w0, h))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_session(sys.argv[1])))
