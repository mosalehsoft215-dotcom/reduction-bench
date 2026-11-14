"""convergence.py: the simulator's quadrature against one with every node count doubled (or more).

One target observation of h1 (companion 291 mas, the hardest case for the channel and time sums) and one of the
first public night, image motion 30 mas rms on every telescope: largest relative change of a squared visibility and
largest change of a closure phase between the settings used for the data (24 time nodes, 48 wavenumber nodes per
channel, 10 x 10 image-motion nodes, coupling as a 4-term polynomial across a channel) and finer or coarser ones."""
import numpy as np

import nights as ng
import skymodel as sm


def set_order(p):
    sm.P_ORDER = p
    sm._CX = np.cos(np.pi * (np.arange(p) + 0.5) / p)
    sm._VINV = np.linalg.inv(sm._CX[:, None] ** np.arange(p)[None, :])


def run(spec, cfg):
    tg = spec["target"]
    m = ng.target_times(spec)[3]
    set_order(cfg[3])
    return ng.observe(ng.system_of(tg, m, True), np.radians(tg["ra"]), np.radians(tg["dec"]), m, 120.0, 1.001,
                      (820.0, 283.0), {t: 30.0 for t in sm.TELS}, tnodes=cfg[0], snodes=cfg[1])


for name in ("h1", "night_2025-07-06"):
    spec = next(s for s in ng.NIGHTS if s["name"] == name)
    fine = run(spec, (48, 96, 14, 6))
    for cfg in ((24, 48, 10, 4), (16, 32, 8, 3)):
        r = run(spec, cfg)
        print(f"{name}: nodes (time, wavenumber, image motion, coupling terms) {cfg} vs (48, 96, 14, 6): "
              f"largest relative change of vis2 {np.max(np.abs(r[3] / fine[3] - 1)):.1e}, "
              f"of closure phase {np.max(np.abs((r[4] - fine[4] + 180) % 360 - 180)):.1e} deg")
