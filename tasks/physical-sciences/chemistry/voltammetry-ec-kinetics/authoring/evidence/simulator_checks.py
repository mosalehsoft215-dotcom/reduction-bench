"""Checks of the disk simulator behind the data (authoring/provenance/microdisk.py).

Steady state (Disk) and transient CVs (transient_cv) at a disk in an
insulating plane, against exact limits and published results for disks. A
separate brute-force solve in cylindrical coordinates is in
microdisk_crosscheck.py.

Run from the bundle root:  python authoring/evidence/simulator_checks.py
Results from the last run are kept in simulator_checks.txt next to this file.
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "provenance"))
import generate_cv_data as g  # noqa: E402
import microdisk as md  # noqa: E402

F = md.FARADAY
f = F / (md.GAS * 298.15)
R = g.RADIUS
AN = dict(e0=g.ANALYTE["e0_vs_ref"], k0=g.ANALYTE["k0"], alpha=g.ANALYTE["alpha"], d_ox=g.ANALYTE["d"],
          d_red=g.ANALYTE["d"], c_ox=g.CONC_TRUE, c_red=0.0, kc=g.ANALYTE["kc"], kb=g.ANALYTE["kb"])


def tri(e1, e2, v, de):
    n_half = int(round(abs(e2 - e1) / de))
    leg = np.linspace(e1, e2, n_half + 1)
    return np.concatenate([leg, leg[-2::-1]]), de / v


def shoup_szabo(tau):
    """Shoup and Szabo (J. Electroanal. Chem. 140 (1982) 237), i / (4 n F D C r) after a
    diffusion-limited step at a disk, tau = 4 D t / r^2; stated accuracy 0.6%."""
    return 0.7854 + 0.8862 / np.sqrt(tau) + 0.2146 * np.exp(-0.7823 / np.sqrt(tau))


def main():
    out = []

    def say(msg):
        print(msg, flush=True)
        out.append(msg)

    t0 = time.time()
    say(f"disk radius {1e4 * R:.1f} um throughout")

    say("1) steady state (Disk): the analyte's parameters without and with R <-> Z")
    dd, cc = AN["d_ox"], AN["c_ox"]
    ilim = -4 * F * dd * cc * R
    disk = md.Disk(R)
    ee = np.linspace(-0.9, -1.5, 121)
    base = {k: AN[k] for k in ("e0", "k0", "alpha", "d_ox", "d_red", "c_ox", "c_red")}
    i_lim = disk.current(np.array([-3.0]), **base)[0]
    say(f"   limiting current / (4 F D C r) = {i_lim / ilim:.10f}")
    i_n = disk.current(ee, **dict(base, k0=1e3, alpha=0.5, e0=-1.19))
    say(f"   Nernstian wave against i_lim / (1 + exp(f (E - E0))): max difference "
        f"{100 * np.abs(i_n - ilim / (1 + np.exp(f * (ee + 1.19)))).max() / abs(ilim):.4f}% of i_lim")
    i_q = disk.current(ee, **base)
    i_o = md.oldham_disk(ee, e0=AN["e0"], k0=AN["k0"], alpha=AN["alpha"], d=dd, c_ox=cc, radius=R)
    say(f"   quasi-reversible against Oldham's approximation (J. Electroanal. Chem. 313 (1991) 3, "
        f"stated accuracy about 0.3%): max difference {100 * np.abs(i_q - i_o).max() / abs(ilim):.3f}% of i_lim")
    i_c = disk.current(ee, **dict(base, kc=AN["kc"], kb=AN["kb"]))
    for nt, ne in ((80, 48), (320, 192)):
        i_f = md.Disk(R, n_theta=nt, n_eta=ne).current(ee, **dict(base, kc=AN["kc"], kb=AN["kb"]))
        say(f"   grid {nt}x{ne} against the 160x96 used for the data (with R <-> Z): max difference "
            f"{100 * np.abs(i_f - i_c).max() / abs(ilim):.4f}% of i_lim")

    say("2) transient CV (transient_cv), reversible couple: forward peak against Aoki et al. "
        "(J. Electroanal. Chem. 171 (1984) 219, a fit to numerical results)")
    for label, d_, win, kw in (("reduction", 0.86e-5, (-0.8, -1.75), dict(e0=-1.19, c_ox=1e-6, c_red=0.0)),
                               ("oxidation", 2.43e-5, (-0.25, 0.45), dict(e0=0.10, c_ox=0.0, c_red=1e-6))):
        for v in (0.05, 0.5, 5.0):
            p = R * np.sqrt(f * v / d_)
            prog, dt = tri(win[0], win[1], v, 0.25e-3)
            i = md.transient_cv(prog, dt, radius=R, k0=100.0, alpha=0.5, d_ox=d_, d_red=d_, **kw)
            half = len(prog) // 2
            ip = (i[:half + 1].max() if kw["c_red"] > 0 else -i[:half + 1].min()) / (4 * F * d_ * 1e-6 * R)
            say(f"   {label}, p = {p:5.2f}: i_p / (4 F D C r) = {ip:.5f}; Aoki {md.aoki_peak(p):.5f} "
                f"({100 * (ip / md.aoki_peak(p) - 1):+.2f}%)")

    say("3) potential step to the diffusion limit against Shoup and Szabo (stated accuracy 0.6%)")
    d_ = 1e-5
    t_end = 5.0 * R * R / d_
    n = 20000
    dt = t_end / n
    prog = np.full(n + 1, -3.0)
    prog[0] = 0.0
    i = md.transient_cv(prog, dt, radius=R, e0=0.0, k0=100.0, alpha=0.5, d_ox=d_, d_red=d_, c_ox=1e-6, c_red=0.0)
    t = dt * np.arange(n + 1)
    tau = 4 * d_ * t / (R * R)
    for tt in (0.05, 0.2, 1.0, 5.0, 20.0):
        k = int(np.argmin(np.abs(tau - tt)))
        val = -i[k] / (4 * F * d_ * 1e-6 * R)
        say(f"   tau = {tau[k]:6.3f}: i / (4 F D C r) = {val:.5f}; Shoup-Szabo {shoup_szabo(tau[k]):.5f} "
            f"({100 * (val / shoup_szabo(tau[k]) - 1):+.2f}%)")

    say("4) a slow sweep of the analyte approaches its steady state (0.5 mV/s, with R <-> Z)")
    prog, dt = tri(-0.8, -1.75, 0.0005, 1e-3)
    i = md.transient_cv(prog, dt, radius=R, **AN)
    iss = disk.current(prog, **AN)
    half = len(prog) // 2
    say(f"   forward sweep against the steady state: max difference "
        f"{100 * np.abs(i[:half] - iss[:half]).max() / abs(ilim):.2f}% of i_lim (the transient excess, "
        f"which falls as the scan slows: 1/sqrt(t) at the plateau)")

    say("5) convergence for the analyte (with R <-> Z and charging), against the settings of the data "
        "(160x96 grid, 0.25 mV steps)")
    for v in (0.05, 10.0):
        ref = None
        for label, de, nt, ne in (("0.25 mV, 160x96", 0.25e-3, 160, 96), ("0.5 mV, 160x96", 0.5e-3, 160, 96),
                                  ("0.25 mV, 240x144", 0.25e-3, 240, 144), ("0.25 mV, 120x64", 0.25e-3, 120, 64)):
            prog, dt = tri(-0.8, -1.75, v, de)
            i = md.transient_cv(prog, dt, radius=R, cdl=g.CDL, n_theta=nt, n_eta=ne, **AN)
            i = i[::int(round(1e-3 / de))]
            if ref is None:
                ref = i
                continue
            say(f"   {v:5g} V/s, {label}: max difference {100 * np.abs(i - ref).max() / np.abs(ref).max():.4f}% "
                f"of the peak")
    say(f"(total {time.time() - t0:.0f} s)")
    with open(os.path.join(os.path.dirname(__file__), "simulator_checks.txt"),
              "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
