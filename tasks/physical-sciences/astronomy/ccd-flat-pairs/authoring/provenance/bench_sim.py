"""Physical model of one flat-pair run on the electro-optical bench.

Charge collection: a flat of mean level E electrons is built up in STEPS slices. In each slice a pixel receives
Poisson(rate * (1 + sum_ij a_ij Q_ij)) electrons, the kernel acting on the charge already stored plus half of the
slice's expected charge (midpoint rule), so that brighter-fatter feeds back as the charge grows. The kernel sums to
zero. A pixel's charge above its own full-well depth spills equally into its two column neighbours. A sensor whose
segment lies at the physical corner of the die has an edge roll-off: the first rows and columns collect less light
by a fraction `edge_rolloff * edge_ratio^k` (no charge moves).

Parallel transfer (row 0 reaches the serial register first): optional deferred charge in the columns, a fraction
`row_eps` of each packet held back and released into the following rows with ratio `row_ratio` per row.
Readout of each row (column 0 first): optional deferred charge in the serial register, a fraction `eps` of each
packet held back and released into the following reads with ratio `ratio` per read; then the output chain, a
static response adu * (1 - k (adu / 60000)^2), the offset, a row offset, white read noise and 16-bit rounding.

Monitor photodiode: charge proportional to the light during the open time, read on the range the electrometer
picks (2, 20, 200 nC full scale); a sensor may carry a calibration factor on one range.
"""
import numpy as np
from numpy.fft import rfft2, irfft2

ROWS, COLS, OVER = 256, 256, 16
STEPS = 64
NEIGHBOURS = [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
RANGES = [(2.0, "2nC"), (20.0, "20nC"), (np.inf, "200nC")]


def area_kernel(a00, weights):
    """Periodic ROWS x COLS kernel: a00 at the origin, the rest shared by the neighbours in proportion."""
    k = np.zeros((ROWS, COLS))
    total = 0.0
    for (dy, dx), w in zip(NEIGHBOURS, weights):
        for sy in ((1,) if dy == 0 else (1, -1)):
            for sx in ((1,) if dx == 0 else (1, -1)):
                k[(sy * dy) % ROWS, (sx * dx) % COLS] += w
                total += w
    k *= -a00 / total
    k[0, 0] = a00
    return k


class Bench:
    def __init__(self, p, rng):
        self.p, self.rng = p, rng
        self.kf = rfft2(area_kernel(p["a00"], p["weights"]))
        self.prnu = 1.0 + rng.normal(0.0, p["prnu"], (ROWS, COLS))
        if p["edge_rolloff"]:          # physical sensor edge at row 0 and column 0: pixels there collect less light
            k = np.arange(ROWS)
            edge = p["edge_rolloff"] * p["edge_ratio"] ** k
            self.prnu *= (1.0 - edge)[:, None] * (1.0 - edge)[None, :]
        self.well = p["full_well"] * (1.0 + rng.normal(0.0, p["well_spread"], (ROWS, COLS)))
        self.columns = rng.normal(0.0, 1.2, COLS + OVER)            # fixed offset pattern, ADU
        y = np.linspace(-1.0, 1.0, ROWS)[:, None]
        x = np.linspace(-1.0, 1.0, COLS)[None, :]
        self.xy = (x, y)

    def collect(self, level, light):
        q = np.zeros((ROWS, COLS))
        rate = level * light / STEPS
        for _ in range(STEPS):
            grow = 1.0 + irfft2(rfft2(q + 0.5 * rate) * self.kf, s=q.shape)
            q = q + self.rng.poisson(np.clip(rate * grow, 0.0, None))
            over = np.clip(q - self.well, 0.0, None)
            if over.any():
                q = q - over + 0.5 * (np.roll(over, 1, 0) + np.roll(over, -1, 0))
        return q

    def read(self, q):
        p, rng = self.p, self.rng
        if p["row_eps"]:
            eps, ratio = p["row_eps"], p["row_ratio"]
            shifted = q * (1.0 - eps)
            held = np.zeros(COLS)
            for i in range(ROWS):
                shifted[i] += held * (1.0 - ratio)
                held = held * ratio + eps * q[i]
            q = shifted
        line = np.concatenate([q, np.zeros((ROWS, OVER))], axis=1)
        if p["defer_eps"]:
            eps, ratio = p["defer_eps"], p["defer_ratio"]
            out = line * (1.0 - eps)
            held = np.zeros(ROWS)
            for j in range(COLS + OVER):
                out[:, j] += held * (1.0 - ratio)
                held = held * ratio + eps * line[:, j]
            line = out
        adu = line / p["gain"]
        adu = adu * (1.0 - p["compress"] * (adu / 60000.0) ** 2)
        adu = adu + p["offset"] + self.columns[None, :] + rng.normal(0.0, 1.0, (ROWS, 1))
        adu = adu + rng.normal(0.0, p["read_e"] / p["gain"], adu.shape)
        return np.clip(np.round(adu), 0, 65535).astype(np.uint16)

    def monitor(self, level, lamp):
        p = self.p
        charge = p["monitor_nc_per_e"] * level * lamp * (1.0 + self.rng.normal(0.0, 2e-4))
        name = next(n for top, n in RANGES if charge < top)
        if name == p.get("miscal_range"):
            charge *= 1.0 + p["miscal"]
        return charge, name

    def run(self):
        p, rng = self.p, self.rng
        n = p["levels"]
        top = p["highest_e"]
        levels = np.concatenate([np.geomspace(p["lowest_e"], 0.08 * top, 8), np.linspace(0.12 * top, top, n - 8)])
        order = rng.permutation(n)
        frames = np.zeros((n, 2, ROWS, COLS + OVER), np.uint16)
        log = []
        clock = 30.0
        x, y = self.xy
        for slot, li in enumerate(order):
            texp = levels[li] / p["lamp_e_per_s"]
            angle = rng.uniform(0.0, 2.0 * np.pi)
            for k in range(2):
                lamp = 1.0 + p["drift"] * np.sin(clock / p["drift_period_s"] + p["drift_phase"]) + rng.normal(0.0, 0.002)
                light = self.prnu * lamp
                if p["pattern_shift"]:
                    side = 0.5 if k == 0 else -0.5
                    light = light * (1.0 + p["pattern_shift"] * side * (np.cos(angle) * x + np.sin(angle) * y))
                frames[li, k] = self.read(self.collect(levels[li], light))
                charge, rname = self.monitor(levels[li], lamp)
                log.append((li, k, round(clock, 1), round(texp, 4), charge, rname))
                clock += texp + 9.0
        zeros = np.stack([self.read(np.zeros((ROWS, COLS))) for _ in range(4)])
        return frames, zeros, log
