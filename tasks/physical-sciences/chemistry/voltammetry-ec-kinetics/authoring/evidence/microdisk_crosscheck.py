"""The microdisk solver against a brute-force solve in cylindrical coordinates.

microdisk.py works in oblate spheroidal coordinates with precomputed Green's
matrices. This script solves the same steady-state problem directly on a
stretched (r, z) grid, refined at the disk edge, with its own finite
volumes, its own treatment of the electrode surface and a far-field
condition for a point sink (c = c_bulk - Q / rho). The two share no code.
Both are run at the analyte's parameters from generate_cv_data.py, with
and without the follow-up chemistry.

Run from the bundle root (numpy and scipy as in environment/Dockerfile):
  python authoring/evidence/microdisk_crosscheck.py
Results from the last run are in microdisk_crosscheck.txt next to this file.
"""

import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "provenance"))
import generate_cv_data as g  # noqa: E402
import microdisk  # noqa: E402

F_RT = microdisk.FARADAY / (microdisk.GAS * g.TEMP)


def faces(s_min, q, top):
    s = [0.0]
    h = s_min
    while s[-1] < top:
        s.append(min(s[-1] + h, top))
        h *= q
    return np.array(s)


class Cylinder:
    """Cell-centred finite volumes on (r, z) in units of the disk radius, with D = 1."""

    def __init__(self, s_min=2e-3, q=1.08, far=60.0):
        inner = 1.0 - faces(s_min, q, 1.0)[::-1]
        outer = 1.0 + faces(s_min, q, far - 1.0)[1:]
        self.rf = np.concatenate([inner, outer])
        self.zf = faces(s_min, q, far)
        self.nr, self.nz = self.rf.size - 1, self.zf.size - 1
        self.rc = 0.5 * (self.rf[:-1] + self.rf[1:])
        self.zc = 0.5 * (self.zf[:-1] + self.zf[1:])
        self.disk = np.where(self.rf[1:] <= 1.0 + 1e-12)[0]      # cells on the electrode
        self.ring = np.pi * (self.rf[1:] ** 2 - self.rf[:-1] ** 2)

    def idx(self, i, j):
        return j * self.nr + i

    def operator(self, lam, far_value_is_zero):
        """Diffusion (+ decay lam) matrix for one species, and the right-hand
        side of the far-field condition for a unit bulk value."""
        nr, nz, rf, zf, rc, zc = self.nr, self.nz, self.rf, self.zf, self.rc, self.zc
        rows, cols, vals = [], [], []
        diag = np.zeros(nr * nz)
        rhs = np.zeros(nr * nz)
        dz = np.diff(zf)
        for j in range(nz):
            for i in range(nr):
                k = self.idx(i, j)
                vol = self.ring[i] * dz[j]
                diag[k] += lam * vol
                # radial faces
                if i + 1 < nr:
                    cond = 2 * np.pi * rf[i + 1] * dz[j] / (rc[i + 1] - rc[i])
                    kk = self.idx(i + 1, j)
                    rows += [k, kk]; cols += [kk, k]; vals += [-cond, -cond]
                    diag[k] += cond; diag[kk] += cond
                else:                                            # outer wall r = far
                    area = 2 * np.pi * rf[-1] * dz[j]
                    rho2 = rf[-1] ** 2 + zc[j] ** 2
                    cond = area * rf[-1] / rho2
                    diag[k] += cond
                    rhs[k] += 0.0 if far_value_is_zero else cond
                # axial faces
                if j + 1 < nz:
                    cond = self.ring[i] / (zc[j + 1] - zc[j])
                    kk = self.idx(i, j + 1)
                    rows += [k, kk]; cols += [kk, k]; vals += [-cond, -cond]
                    diag[k] += cond; diag[kk] += cond
                else:                                            # top wall z = far
                    rho2 = rc[i] ** 2 + zf[-1] ** 2
                    cond = self.ring[i] * zf[-1] / rho2
                    diag[k] += cond
                    rhs[k] += 0.0 if far_value_is_zero else cond
        rows += list(range(nr * nz)); cols += list(range(nr * nz)); vals += list(diag)
        return coo_matrix((vals, (rows, cols)), shape=(nr * nz, nr * nz)).tocsr(), rhs

    def current(self, e, *, e0, k0, alpha, d, c_ox, radius, kc=0.0, kb=0.0):
        """Ratio of the current to 4 F D C a, for a reduction with no R in the bulk."""
        a = radius
        n = self.nr * self.nz
        nd = self.disk.size
        kap0 = k0 * a / d
        lam = (kc + kb) * a * a / d
        ao, bo = self.operator(0.0, False)          # O, bulk 1
        as_, _ = self.operator(0.0, True)           # S = R + Z, bulk 0
        aw, _ = self.operator(lam, True)            # W' = (kc R - kb Z)/(kc + kb), bulk 0
        th_b = kb / (kc + kb) if kc > 0 else 1.0     # R = th_b S + W'
        th_f = 1.0 - th_b
        dz0 = self.zf[1] - self.zf[0]
        out = []
        for ee in np.atleast_1d(e):
            eta = F_RT * (ee - e0)
            e1, e2 = kap0 * np.exp(-alpha * eta), kap0 * np.exp((1 - alpha) * eta)
            # unknowns: O, S, W' in the cells, then O_s, S_s, W'_s at the disk cells
            rows, cols, vals = [], [], []
            rhs = np.zeros(3 * n + 3 * nd)
            for blk, mat in enumerate((ao, as_, aw)):
                m = mat.tocoo()
                rows += list(m.row + blk * n); cols += list(m.col + blk * n); vals += list(m.data)
            rhs[:n] = bo
            for q, i in enumerate(self.disk):
                k = self.idx(i, 0)
                area = self.ring[i]
                cs = area / (0.5 * dz0)                   # cell centre to surface
                so, ss, sw = 3 * n + q, 3 * n + nd + q, 3 * n + 2 * nd + q
                for blk, sidx in ((0, so), (1, ss), (2, sw)):
                    # cell: flux from the surface node into the cell
                    rows += [blk * n + k, blk * n + k]; cols += [blk * n + k, sidx]; vals += [cs, -cs]
                # surface nodes: cs (c_cell - c_s) = flux into the electrode (O), out of it (S, W')
                # J = area [e1 O_s - e2 (th_b S_s + W'_s)]
                jo = [(so, area * e1), (ss, -area * e2 * th_b), (sw, -area * e2)]
                # O:  cs (O_cell - O_s) - J = 0
                rows += [so, so]; cols += [k, so]; vals += [cs, -cs]
                for c_, v_ in jo:
                    rows.append(so); cols.append(c_); vals.append(-v_)
                # S:  cs (S_cell - S_s) + J = 0
                rows += [ss, ss]; cols += [n + k, ss]; vals += [cs, -cs]
                for c_, v_ in jo:
                    rows.append(ss); cols.append(c_); vals.append(v_)
                # W': cs (W_cell - W_s) + th_f J = 0
                rows += [sw, sw]; cols += [2 * n + k, sw]; vals += [cs, -cs]
                for c_, v_ in jo:
                    rows.append(sw); cols.append(c_); vals.append(th_f * v_)
            mat = coo_matrix((vals, (rows, cols)), shape=(3 * n + 3 * nd,) * 2).tocsc()
            x = spsolve(mat, rhs)
            o_s, s_s, w_s = x[3 * n:3 * n + nd], x[3 * n + nd:3 * n + 2 * nd], x[3 * n + 2 * nd:]
            flux = self.ring[self.disk] * (e1 * o_s - e2 * (th_b * s_s + w_s))
            out.append(flux.sum() / 4.0)
        return np.array(out)


def main():
    lines = []

    def say(msg):
        print(msg, flush=True)
        lines.append(msg)

    p = g.ANALYTE
    a, d, c = g.RADIUS, p["d"], g.CONC_TRUE
    ilim = -4 * microdisk.FARADAY * d * c * a
    e = np.array([-1.05, -1.12, -1.17, -1.20, -1.23, -1.27, -1.33, -1.45, -1.75])
    disk = microdisk.Disk(a)
    for label, kc, kb in (("no chemistry", 0.0, 0.0),
                          (f"R <-> Z, kc {p['kc']} 1/s, kb {p['kb']} 1/s", p["kc"], p["kb"])):
        spher = disk.current(e, e0=p["e0_vs_ref"], k0=p["k0"], alpha=p["alpha"], d_ox=d,
                             d_red=d, c_ox=c, c_red=0.0, kc=kc, kb=kb) / ilim
        say(f"{label}: current / (4 F D C a) at the analyte's k0 and alpha")
        prev = None
        for s_min, q in ((4e-3, 1.12), (2e-3, 1.08), (1e-3, 1.06)):
            t0 = time.time()
            cyl = Cylinder(s_min, q).current(e, e0=p["e0_vs_ref"], k0=p["k0"], alpha=p["alpha"],
                                             d=d, c_ox=c, radius=a, kc=kc, kb=kb)
            diff = cyl - spher
            say(f"  (r, z) grid {Cylinder(s_min, q).nr}x{Cylinder(s_min, q).nz}: max |difference| "
                f"{100 * np.abs(diff).max():.3f}% of the limiting current ({time.time() - t0:.0f} s)")
            prev = cyl
        for ee, s_, c_ in zip(e, spher, prev):
            say(f"    E {ee:+.3f} V: spheroidal {s_:.5f}, cylindrical {c_:.5f}")
    text = "\n".join(lines)
    with open(os.path.join(HERE, "microdisk_crosscheck.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
