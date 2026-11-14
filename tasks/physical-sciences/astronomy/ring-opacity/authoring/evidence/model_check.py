"""The two parts of the light at each crossing, and the reference's shortcuts against the generator.

    PYTHONPATH=../../solution:../provenance python model_check.py   (in the task image)  -> model_check.txt
Part 1: at the middle of each of the sixteen crossings, the direct beam and the diffracted light of the
generator (noise-free, fractions of the star's flux), and the share of the diffracted light that was
diffracted more than once.
Part 2: the reference computes the direct beam with 128 nodes per exposure and only near the edges, and the
diffracted light from fewer ring elements and orders than the generator. Both are compared with the
generator's noise-free curves at the generator's ring, in the windows that the reference fits.
"""
from pathlib import Path

import numpy as np

import generate_ring as G
import ring_fit as R
from ring_model import halo, near_edges, patterns, ring_plane, star_flux

HERE = Path(__file__).resolve().parent
r = G.ring()
table = patterns(R.R_MIN_KM, R.R_MAX_KM, R.SLOPE, G.LAM_KM, G.LAM_W, G.DIST_KM)
lines = ["station, time of the crossing (s), ring longitude (deg), radial width (km), direct beam, diffracted light, "
         "share diffracted more than once"]
far = near = light_diff = 0.0
for st, t, beam, light in G.clean():
    mid = t + st["exposure"] / 2
    rho, th = ring_plane(st["xi0"] + st["vxi"] * mid, st["eta0"] + st["veta"] * mid, r["centre"], r["sin_b"], r["pa"])
    win = np.zeros(t.size, bool)
    for i in np.flatnonzero(np.diff(np.sign(rho - r["a"]))):
        x = r["tau_los"](th[i]) / 2
        once = x * np.exp(-2 * x) / (np.exp(-x) - np.exp(-2 * x))
        lines.append(f"{st['name']}  {mid[i]:7.2f}  {np.degrees(th[i]) % 360:6.1f}  {G.W0 * (1 - G.Q_WIDTH * np.cos(th[i] - G.THETA_P)):5.2f}  "
                     f"{np.exp(-2 * x):.4f}  {light[i]:.4f}  {1 - once:.2f}")
        win |= np.abs(mid - mid[i]) < R.HALF
    rows = near_edges(t[win], st, r, R.NEAR_KM)
    d = np.abs(star_flux(t[win], st, r, G.LAM_KM, G.LAM_W, G.DIST_KM, rows=rows) - beam[win])
    far, near = max(far, d[~rows].max()), max(near, d[rows].max())
    mine = halo(st["xi0"] + st["vxi"] * mid[win], st["eta0"] + st["veta"] * mid[win], r, table, orders=R.ORDERS)
    light_diff = max(light_diff, np.abs(mine - light[win]).max())
lines += ["", "reference against the generator in the fit windows, largest |difference| in units of the star's flux:",
          f"direct beam, exposures within {R.NEAR_KM:g} km of an edge (diffraction integral, 128 nodes): {near:.1e}",
          f"direct beam, the other exposures (geometric optics): {far:.1e}",
          f"diffracted light (1440 x 6 ring elements, {R.ORDERS} orders, against 2160 x 8 and 12): {light_diff:.1e}"]
(HERE / "model_check.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
print("\n".join(lines))
