"""Reference analysis: the electron source from the plasma-lens focus scan, shot by shot.

1. Spectrometer: the background (median of the outer pixels) is removed from every profile, counts
   per pixel become counts per MeV with the pixel widths of the energy axis, each shot is normalised
   to unit area and the shots are averaged. This gives the mean kinetic energy, the measured rms
   relative spread of the run and, from the shot means, the shot-to-shot jitter of the energy. The
   spectrometer has no focusing, so the measured spread contains the beam's own horizontal spot on
   its screen (R11 x + R12 x' from the lens entrance, plus the 0.10 mm screen resolution), converted
   to energy with the local dispersion of the energy axis. That blur depends on the source the scan
   measures, so it is recomputed from the fitted source in every round of step 3.
   The momentum distribution within a shot is not Gaussian. Its shape comes from windowed central
   moments of every shot in screen pixels: the spot is a Gaussian in the screen coordinate, so it
   adds to the second cumulant and leaves the third unchanged. The shot's own screen distribution
   (spread and skewness, skew-normal) is mapped through the energy axis into momentum, and its
   nodes replace Gaussian nodes in step 3.
2. Scan profiles: every shot gets a straight background from its outer pixels and is shifted exactly
   (Fourier shift) onto its centroid. The shots are not averaged.
3. Fit. Two things show in the shots and are part of the model.
   - Shots at the same set current have different widths, far more than the energy jitter of the
     spectrometer explains: the focusing strength itself varies from shot to shot. Every shot
     therefore gets its own offset of the ratio momentum / lens current, found by a one-dimensional
     search for that shot. The offsets are kept free of a mean and of a trend with the set current,
     which would stand for another source distance and another divergence.
   - The profiles are not symmetric: before the focus they lean to one side, beyond it to the
     other. A beam that enters the lens off its axis is steered, and each momentum by its own
     amount, so the momentum slices of a shot land side by side. The model displaces every slice
     along the centroid ray M11 x of an entry offset x (an entry angle acts in nearly the same way
     over this scan and is carried by x).
   - The beam's pointing changes from shot to shot (setup.txt). A ray leaving the source at angle th crosses
     the lens entrance at (d th, th), so every shot is steered by its own amount on top of the lens offset:
     each shot gets its own pointing next to its focusing offset. Near focus a shot can be matched with
     its focus on either side and the steering reversed, so there the two are searched together on a grid
     before both are refined. The rms pointing follows from the shots' pointings and their fit errors by
     maximum likelihood. On the spectrometer screen a shot's pointing moves its spectrum by
     (R11 d + R12) th, so the scatter of the spectrometer's shot means holds the pointing as well as the
     energy jitter; the measured rms pointing is removed from the run's spread in quadrature.
   - The lens is a discharge in a plasma of about 1e17 cm^-3 (setup.txt), and the shots differ in charge.
     The beam's own wake in that plasma focuses it in proportion to its charge, about its own centroid: it
     adds kp q / q_mean to the focusing of the envelope of every momentum slice and does not steer the
     centroid rays. The profiles alone do not show it, because the shot offsets absorb it, but then the
     offsets follow the shot's counts above background (its charge) and their mean, which the jitter of the
     discharge does not have, is the wake's focusing at the mean charge. In every round the offsets are
     regressed on the set-current trend and on q / k, and the charge coefficient goes into kp, with no
     focusing at zero charge. kp is reported in T/m.
   Each shot is the mixture over its momentum distribution (shape and rms from step 1, without the energy
   jitter) of Gaussians with the slice's own size and place, centred like the data and blurred by
   the screen (step 4); amplitude and offset of every shot are linear and solved exactly. Free for
   all shots together: rms source size, rms divergence, source distance and x. The shot offsets and
   these four are found in turn until they stop changing.
4. Screen. The crystal is turned by 45 degrees and the camera looks along its normal, so pixel
   calibration and knife edge are lengths on the crystal face; in the beam's horizontal coordinate
   they are shorter by cos 45. The knife edge gives the optics (Gaussian, rms = 10-90 % rise /
   2.5631). The crystal is 50 um thick and an electron emits along its whole path: seen along the
   normal each electron is a line of thickness x tan 45 on the face, thickness x sin 45 in the beam's
   coordinate. Every profile is therefore also convolved with a box of that width.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, least_squares, nnls

from lens_model import C_GEV, ME, centre, hermite, m11_m12, momentum_gev, skew_nodes

G_PER_A = 0.2000                                                  # setup.txt
FACE_PIX, EDGE_RISE = 6.80e-6, 43.5e-6                            # setup.txt, both measured on the crystal face
CRYSTAL, TILT = 50.0e-6, np.radians(45.0)                         # setup.txt
# The camera looks along the crystal normal and the crystal is turned by TILT about the vertical axis, so
# a horizontal length u on the face is u cos(TILT) in the beam's horizontal coordinate.
PIX = FACE_PIX * np.cos(TILT)
RES = EDGE_RISE * np.cos(TILT) / 2.5631    # optics: the 10-90 % rise of a Gaussian edge is 2.5631 sigma
# Depth of the light source. An electron emits along its whole path in the crystal. Seen along the normal
# that path covers CRYSTAL tan(TILT) on the face (light leaving along the normal is not refracted), which
# is CRYSTAL sin(TILT) in the beam's coordinate. Every electron is therefore imaged as a line of that
# length: a box that widens every profile.
BOX = CRYSTAL * np.tan(TILT) * np.cos(TILT)
GL_X, GL_W = np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)]) / 2, np.array([5.0, 8.0, 5.0]) / 18   # box, 3-point Gauss-Legendre
R11_SPEC, R12_SPEC, SPEC_RES_MM, SPEC_MM_PER_PX = 0.92, 1.85, 0.10, 0.10
KEYS = ("source_size_um", "divergence_mrad", "source_distance_mm", "energy_MeV", "energy_spread_percent")
XQ, WQ = hermite(21)                                              # momentum nodes within a shot (updated to the measured shape)
PASSIVE_MODE = os.environ.get("LENS_FIT_PASSIVE_MODE", "charge")  # authoring controls only
COUNTS = []                                                       # counts of every scan shot: its charge, up to one factor


def load_csv(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)


def spectrum(d):
    raw = load_csv(d / "spectrometer.csv")[:, 1:]
    cal = load_csv(d / "energy_axis.csv")[:, 1]
    edges = np.r_[cal[0] - 0.5 * (cal[1] - cal[0]), 0.5 * (cal[1:] + cal[:-1]), cal[-1] + 0.5 * (cal[-1] - cal[-2])]
    width = np.diff(edges)
    dens, means, m2s, m3s, cs = [], [], [], [], []
    px = np.arange(raw.shape[1], dtype=float)
    for y in raw:
        b = np.median(np.r_[y[:40], y[-40:]])
        f = (y - b) / width
        f = f / np.sum(f * width)
        dens.append(f)
        means.append(np.sum(f * width * cal))
        z = y - b                                   # counts per pixel: the screen coordinate carries the blur additively
        c = np.sum(z * px) / np.sum(z)
        sd = np.sqrt(np.sum(z * (px - c) ** 2) / np.sum(z))
        for _ in range(3):
            win = np.abs(px - c) < 6 * sd
            c = np.sum(z[win] * px[win]) / np.sum(z[win])
            sd = np.sqrt(np.sum(z[win] * (px[win] - c) ** 2) / np.sum(z[win]))
        m2s.append(sd ** 2)
        m3s.append(np.sum(z[win] * (px[win] - c) ** 3) / np.sum(z[win]))
        cs.append(c)
    global SPEC_M2, SPEC_M3, SPEC_C, SPEC_CAL
    SPEC_M2, SPEC_M3, SPEC_C, SPEC_CAL = float(np.mean(m2s)), float(np.mean(m3s)), float(np.mean(cs)), cal
    f = np.mean(dens, axis=0)
    w = f * width
    mean = np.sum(w * cal) / np.sum(w)
    rel = np.sqrt(np.sum(w * (cal - mean) ** 2) / np.sum(w)) / mean     # measured, blur included
    jit = np.std(means, ddof=1) / mean                                  # shot-to-shot jitter of the mean
    e_tot = mean + ME
    kin = (e_tot ** 2 - ME ** 2) / (e_tot * mean)                    # dT/T = kin * dp/p
    j = np.interp(mean, cal, np.arange(cal.size))
    disp = np.interp(j, np.arange(cal.size - 1) + 0.5, np.diff(np.log(cal))) / SPEC_MM_PER_PX   # dlnT/ds, per mm
    return mean, rel, kin, disp, jit


def blur_rel(size, div, dist, disp):
    """relative kinetic-energy blur of the spectrometer: the beam's horizontal spot on its screen (m -> mm)"""
    s2 = (R11_SPEC ** 2 * size ** 2 + (R11_SPEC * dist + R12_SPEC) ** 2 * div ** 2) * 1e6
    return np.sqrt(s2 + SPEC_RES_MM ** 2) * disp


def shots(d):
    """every shot: set current, positions (m, from its centroid), counts above background, weights"""
    a = load_csv(d / "focus_scan.csv")
    cur, y = a[:, 0], a[:, 2:]
    n = y.shape[1]
    px = np.arange(n)
    xc = (px - n / 2 + 0.5) * PIX
    outer = np.r_[0:150, n - 150:n]
    freq = np.fft.rfftfreq(n)
    out = []
    for I, yy in zip(cur, y):
        c = np.polyfit(px[outer], yy[outer], 1)
        z = yy - np.polyval(c, px)
        read = np.std(z[outer])
        core = np.abs(px - np.argmax(np.convolve(z, np.ones(9) / 9, "same"))) < 300
        x0 = np.sum(z[core] * xc[core]) / np.sum(z[core])
        rms = np.sqrt(max(np.sum(z[core] * (xc[core] - x0) ** 2) / np.sum(z[core]), (3 * PIX) ** 2))
        zs = np.fft.irfft(np.fft.rfft(z) * np.exp(2j * np.pi * freq * x0 / PIX), n)   # centred, no smoothing
        win = np.abs(xc) < np.clip(6 * rms, 60 * PIX, 420 * PIX)
        var = read ** 2 + np.clip(np.convolve(zs, np.ones(5) / 5, "same"), 0, None)   # read noise + counting
        out.append((I, xc[win], zs[win], 1.0 / np.sqrt(var[win])))
        COUNTS.append(float(np.sum(z[core])))
    return out


def shot_shapes(x, k_nom, dl, q, within, angle=False, th=0.0, kp=0.0):
    """centred unit-area profiles of one shot for trial offsets dl of its momentum / current ratio;
    q = (size um, divergence mrad, distance mm, entry offset um); returns (len(dl), len(x)).
    th is the shot's pointing at the source in urad (scalar or one value per row of dl): it moves the
    centroid at the lens entrance by dist*th and th. angle=True reads q[3] as an entry angle in urad
    instead (an authoring variant)"""
    size, div, dist, off = q[0] * 1e-6, q[1] * 1e-3, q[2] * 1e-3, q[3] * 1e-6
    kc = k_nom / (1 + dl[:, None] + XQ[None, :] * within)
    # the beam's own wake in the discharge plasma adds kp; it focuses about the beam's centroid, so it enters
    # the sizes and not the steering of the centroid rays
    m11, m12 = m11_m12(kc + kp / (1 + XQ[None, :] * within))
    s2 = m11 ** 2 * size ** 2 + (m11 * dist + m12) ** 2 * div ** 2 + RES ** 2
    m11, m12 = m11_m12(kc)
    thr = np.broadcast_to(np.asarray(th, float) * 1e-6, dl.shape)[:, None]
    c = (m12 if angle else m11) * off + (m11 * dist + m12) * thr
    c = c - np.sum(WQ * c, axis=1, keepdims=True)
    out = 0.0
    for node, wt in zip(GL_X * BOX, GL_W):                       # the box of the crystal depth
        out = out + wt * np.sum(WQ[None, :, None] * np.exp(-(x[None, None, :] - c[:, :, None] - node) ** 2 / (2 * s2[:, :, None]))
                                / np.sqrt(2 * np.pi * s2[:, :, None]), axis=1)
    return out


def linear(z, w, g):
    """best amplitude and offset of every row of g (model shapes) against the shot, and its chi2"""
    w2 = w * w
    sgg, sg1, s11 = (g * g) @ w2, g @ w2, w2.sum()
    sgz, s1z, szz = g @ (w2 * z), np.sum(w2 * z), np.sum(w2 * z * z)
    det = sgg * s11 - sg1 ** 2
    amp, off = (sgz * s11 - sg1 * s1z) / det, (sgg * s1z - sg1 * sgz) / det
    return szz - amp * sgz - off * s1z, amp, off


def shot_offset(x, z, w, k_nom, q, within, start, wide, angle=False, th=0.0, kp=0.0):
    """offset of the shot's momentum / current ratio: grid search, zoomed twice, then a parabola"""
    d0 = start
    for half, npt in (((0.024, 25),) if wide else ()) + ((0.002, 9), (0.0005, 9)):
        dl = d0 + np.linspace(-half, half, npt)
        chi = linear(z, w, shot_shapes(x, k_nom, dl, q, within, angle, th=th, kp=kp))[0]
        i = int(np.argmin(chi))
        d0 = dl[i]
    if 0 < i < npt - 1 and chi[i - 1] - 2 * chi[i] + chi[i + 1] > 0:
        d0 += 0.5 * (dl[1] - dl[0]) * (chi[i - 1] - chi[i + 1]) / (chi[i - 1] - 2 * chi[i] + chi[i + 1])
    return d0


def shot_pointing(x, z, w, k_nom, q, within, dl, start, wide, kp=0.0):
    """pointing of the shot (urad) and its standard error: grid search, zoomed twice, then a parabola"""
    t0 = start
    for half, npt in (((2400.0, 25),) if wide else ()) + ((240.0, 9), (60.0, 9)):
        th = t0 + np.linspace(-half, half, npt)
        chi = linear(z, w, shot_shapes(x, k_nom, np.full(npt, dl), q, within, th=th, kp=kp))[0]
        i = int(np.argmin(chi))
        t0 = th[i]
    h = th[1] - th[0]
    curv = (chi[max(i - 1, 0)] - 2 * chi[i] + chi[min(i + 1, npt - 1)]) / h ** 2 if 0 < i < npt - 1 else 0.0
    if curv > 0:
        t0 += 0.5 * h * (chi[i - 1] - chi[i + 1]) / (chi[i - 1] - 2 * chi[i] + chi[i + 1])
    se = np.sqrt(2.0 / curv) if curv > 0 else np.inf          # chi2 = min + (th - t0)^2 / se^2
    return t0, se


def shot_grid(x, z, w, k_nom, q, within, dl, kp=0.0):
    """near focus a shot can be matched with its focus on either side and the steering reversed, so the
    focusing offset and the pointing are searched together on a grid before they are refined"""
    dls, ths = np.meshgrid(dl + np.linspace(-0.02, 0.02, 21), np.linspace(-2400.0, 2400.0, 25), indexing="ij")
    chi = linear(z, w, shot_shapes(x, k_nom, dls.ravel(), q, within, th=ths.ravel(), kp=kp))[0]
    i = int(np.argmin(chi))
    return dls.ravel()[i], ths.ravel()[i]


def shot_joint(x, z, w, k_nom, q, within, dl, th, kp=0.0):
    """joint refinement of the shot's focusing offset and pointing (amplitude and offset solved exactly);
    returns (dl, th, standard error of th)"""
    def res(v):
        g = shot_shapes(x, k_nom, np.array([v[0]]), q, within, th=v[1], kp=kp)
        _, amp, off = linear(z, w, g)
        return w * (z - amp[0] * g[0] - off[0])
    f = least_squares(res, [dl, th], x_scale=[1e-3, 100.0], diff_step=[1e-6, 1e-4], max_nfev=60)
    try:
        cov = np.linalg.inv(f.jac.T @ f.jac)
        se = float(np.sqrt(cov[1, 1]))
    except np.linalg.LinAlgError:
        se = np.inf
    return f.x[0], f.x[1], se


def pointing_rms(th, se, scale=1.0):
    """rms of the pointing over the shots, by maximum likelihood with each shot's own fit error"""
    th = th - th.mean()
    v = (se * scale) ** 2
    from scipy.optimize import minimize_scalar
    f = lambda s: np.sum(np.log(s * s + v) + th ** 2 / (s * s + v))
    return minimize_scalar(f, bounds=(1.0, 5000.0), method="bounded").x


def residuals(q, data, k_nom, delta, within, angle=False, point=None, kps=None):
    out = []
    point = np.zeros(len(data)) if point is None else point
    kps = np.zeros(len(data)) if kps is None else kps
    for (I, x, z, w), k, dl, th, kp in zip(data, k_nom, delta, point, kps):
        g = shot_shapes(x, k, np.array([dl]), q, within, angle, th=th, kp=kp)
        _, amp, off = linear(z, w, g)
        out.append(w * (z - amp[0] * g[0] - off[0]))
    return np.concatenate(out)


def start_values(data, p_gev):
    """starting point (um, mrad, mm): Gaussian widths of the shots with single-energy optics. For a trial
    source distance the squared size and divergence are linear, so the distance is scanned."""
    cur = np.array([s[0] for s in data])
    v = []
    for I, x, z, w in data:
        g = lambda t, a, s, b: a * np.exp(-t ** 2 / (2 * s ** 2)) + b
        q, _ = curve_fit(g, x, z, p0=[z.max(), 50e-6, 0.0])
        v.append(max(q[1] ** 2 - RES ** 2 - BOX ** 2 / 12, 1e-12))
    v = np.array(v)
    currents = np.unique(cur)
    y = np.array([np.mean(v[cur == I]) for I in currents])
    m11, m12 = m11_m12(C_GEV * G_PER_A * currents / p_gev)
    best = None
    for dist in np.linspace(0.010, 0.200, 381):
        A = np.vstack([m11 ** 2, (m11 * dist + m12) ** 2]).T / y[:, None]
        c, rn = nnls(A, np.ones_like(y))
        if best is None or rn < best[0]:
            best = (rn, c, dist)
    _, c, dist = best
    return max(np.sqrt(c[0]) * 1e6, 0.3), max(np.sqrt(c[1]) * 1e3, 0.1), dist * 1e3


SKEW = [0.0]


def within_shape(b_px):
    """standardised momentum distribution within a shot. On the spectrometer screen the blur is a Gaussian in
    the screen coordinate, so it adds to the second cumulant there and leaves the third unchanged: the
    shot's own spread and skewness on the screen follow from the measured cumulants. That screen
    distribution is mapped through the (slightly nonlinear) energy axis into momentum."""
    var = max(SPEC_M2 - b_px ** 2, 1.0)
    g1 = float(np.clip(SPEC_M3 / var ** 1.5, -0.98, 0.98))
    u, w = skew_nodes(g1, n=49, span=6.0)
    t = np.interp(SPEC_C + np.sqrt(var) * u, np.arange(SPEC_CAL.size), SPEC_CAL)
    pm = momentum_gev(t)
    d = pm / np.sum(w * pm) - 1
    sd = np.sqrt(np.sum(w * d * d))
    SKEW[0] = float(np.sum(w * d ** 3) / sd ** 3)
    return d / sd, w


def solve(d, steering=True, per_shot=True, rounds=16, angle=False, free_within=False, remove_blur=True, energy=None, shape=True, pointing=True):
    """full analysis. The keyword variants are authoring controls: steering=False keeps
    the beam on the lens axis; per_shot=False gives every shot its set current; angle=True steers by an
    entry angle instead of an offset; free_within=True fits the within-shot spread from the profiles;
    remove_blur=False takes the spectrometer spread as measured; energy overrides the mean energy."""
    t_mean, measured_t, kin, disp, jit_t = spectrum(d)
    if energy is not None:
        t_mean = energy
    p = momentum_gev(t_mean)
    COUNTS.clear()
    data = shots(d)
    k_nom = C_GEV * G_PER_A * np.array([s[0] for s in data]) / p
    # the shot's charge, relative to the run mean (counts above background are proportional to it)
    qrel = np.array(COUNTS) / np.mean(COUNTS)
    if PASSIVE_MODE == "constant":                # authoring control: a charge-blind extra focusing
        qrel = np.ones_like(qrel)
    elif PASSIVE_MODE == "centred":               # control: only the charge's deviation from the mean focuses
        qrel = qrel - 1.0
    passive = not os.environ.get("LENS_FIT_NO_PASSIVE")
    kpass = 0.0                      # extra focusing of the beam's own wake at the mean charge, k units
    jit_e = jit_t / kin
    start = start_values(data, p)
    q = np.array([0.6 * start[0], start[1], start[2], (500.0 if angle else 100.0) if steering else 0.0])
    nfree = 4 if steering else 3
    lo, hi = [0.05, 0.05, 1.0, -5000.0][:nfree], [100.0, 20.0, 500.0, 5000.0][:nfree]
    scale = [1.0, 0.1, 1.0, 300.0 if angle else 50.0][:nfree]
    if free_within:
        lo, hi, scale = lo + [1e-4], hi + [0.05], scale + [0.002]
    delta = np.zeros(len(data))
    point, point_se = np.zeros(len(data)), np.full(len(data), np.inf)
    trend = np.vstack([np.ones(len(data)), np.log(k_nom) - np.log(k_nom).mean()]).T
    # a shot's focusing offset dl moves k_nom / (1 + dl); the wake adds kpass * qrel to k, so an error in kpass
    # leaves dl ~ -(error) qrel / k_nom: the charge column of the regression measures it
    trendq = np.vstack([trend.T, qrel / k_nom]).T

    def detrend(delta, kpass):
        if not passive:
            return delta - trend @ np.linalg.lstsq(trend, delta, rcond=None)[0], kpass
        coef = np.linalg.lstsq(trendq, delta, rcond=None)[0]
        return delta - trendq @ coef, kpass - coef[2]
    spread_t = measured_t
    within = np.sqrt(max((spread_t / kin) ** 2 - jit_e ** 2, 1e-10))
    kprev = np.inf
    for it in range(rounds):
        if not free_within:
            within = np.sqrt(max((spread_t / kin) ** 2 - jit_e ** 2, 1e-10))
        if per_shot:
            delta = np.array([shot_offset(x, z, w, k, q, within, dl, it == 0, angle, th, kpass * qq)
                              for (I, x, z, w), k, dl, th, qq in zip(data, k_nom, delta, point, qrel)])
            # the jitter is random about the set point: no mean and no trend with the set current (a mean
            # would stand for another source distance, a trend for another divergence); the part that follows
            # the shot's charge is the wake's focusing and goes into kpass
            delta, kpass = detrend(delta, kpass)
        if pointing and it >= 1:
            # the beam's pointing at the source changes from shot to shot and steers each shot by its own amount;
            # its mean is part of the lens offset q[3]. A wide scan first, then focusing offset and pointing jointly.
            starts = []
            for (I, x, z, w), k, dl, th, qq in zip(data, k_nom, delta, point, qrel):
                if x.size <= 400 and not os.environ.get("LENS_FIT_NO_GRID"):   # near focus: both mirror solutions
                    starts.append(shot_grid(x, z, w, k, q, within, dl, kpass * qq))
                elif it <= 2:
                    starts.append((dl, shot_pointing(x, z, w, k, q, within, dl, th, True, kpass * qq)[0]))
                else:
                    starts.append((dl, th))
            ps = [shot_joint(x, z, w, k, q, within, a, b, kpass * qq)
                  for (I, x, z, w), k, (a, b), qq in zip(data, k_nom, starts, qrel)]
            delta = np.array([a for a, b, c in ps])
            point, point_se = np.array([b for a, b, c in ps]), np.array([c for a, b, c in ps])
            delta, kpass = detrend(delta, kpass)
            point -= point.mean()
        if free_within:
            f = least_squares(lambda v: residuals(np.r_[v[:nfree], q[nfree:]], data, k_nom, delta, v[-1], angle, point, kpass * qrel),
                              np.r_[q[:nfree], within], bounds=(lo, hi), x_scale=scale, max_nfev=40)
            within = f.x[-1]
        else:
            f = least_squares(lambda v: residuals(np.r_[v, q[nfree:]], data, k_nom, delta, within, angle, point, kpass * qrel), q[:nfree],
                              bounds=(lo, hi), x_scale=scale, max_nfev=40)
        step = np.max(np.abs(f.x[:3] / q[:3] - 1))
        q = np.r_[f.x[:nfree], q[nfree:]]
        blur = blur_rel(q[0] * 1e-6, q[1] * 1e-3, q[2] * 1e-3, disp)
        if shape:
            b_px = np.sqrt((R11_SPEC ** 2 * (q[0] * 1e-6) ** 2 + (R11_SPEC * q[2] * 1e-3 + R12_SPEC) ** 2 * (q[1] * 1e-3) ** 2) * 1e6
                           + SPEC_RES_MM ** 2) / SPEC_MM_PER_PX
            global XQ, WQ
            XQ, WQ = within_shape(b_px)
        if remove_blur:
            spread_t = np.sqrt(max(measured_t ** 2 - blur ** 2, 1e-8))
        if step < 2e-4 and it > 2 and abs(kpass - kprev) < 0.005:
            break
        kprev = kpass
    if free_within:
        spread_t = np.hypot(within, jit_e) * kin
    chi2 = 2 * f.cost / (f.fun.size - f.x.size - (4 if pointing else 3) * len(data))
    # spectrometer: a shot's pointing moves its spectrum by (R11 d + R12) th on the screen, so the scatter of the
    # shot means holds the pointing as well as the energy jitter. The scan measures the pointing, which is removed.
    if os.environ.get("LENS_FIT_DUMP"):
        np.savez(os.environ["LENS_FIT_DUMP"], delta=delta, point=point, point_se=point_se, q=q, qrel=qrel, kpass=kpass)
    sig_th = pointing_rms(point, point_se, np.sqrt(max(chi2, 1.0))) if pointing else 0.0
    point_rel = (R11_SPEC * q[2] * 1e-3 + R12_SPEC) * sig_th * 1e-6 * 1e3 * disp
    spread_run = np.sqrt(max(spread_t ** 2 - point_rel ** 2, 1e-8))
    return dict(q=q, chi2=chi2, t_mean=t_mean, measured_t=measured_t, spread_t=spread_run, spread_with_pointing=spread_t,
                kin=kin, blur=blur, p=p, jit_e=jit_e, within=within, delta=delta, rounds=it + 1, start=start, data=data,
                k_nom=k_nom, point=point, point_se=point_se, sig_th=sig_th, point_rel=point_rel,
                passive_T_per_m=kpass * p / C_GEV)


def main(data_dir, out_path):
    r = solve(Path(data_dir), shape=not os.environ.get("LENS_FIT_GAUSSIAN"),      # variant switches for authoring controls
              pointing=not os.environ.get("LENS_FIT_NO_POINTING"))
    q, t_mean, spread_t = r["q"], r["t_mean"], r["spread_t"]
    out = dict(zip(KEYS, (round(float(q[0]), 5), round(float(q[1]), 5), round(float(q[2]), 4), round(float(t_mean), 3),
                          round(100 * float(spread_t), 4))))
    k_mid = np.median(r["k_nom"])
    entry = (q[3] * 1e-6, 0.0)
    disp_mid = (centre(k_mid / 1.001, entry) - centre(k_mid / 0.999, entry)) / 0.002     # d(centroid)/d(delta)
    gam = 1.0 + t_mean / ME
    out.update(normalised_emittance_mm_mrad=round(float(q[0] * q[1] * 1e-3 * np.sqrt(gam ** 2 - 1)), 5),
               measured_spectrometer_spread_percent=round(100 * float(r["measured_t"]), 4),
               spectrometer_blur_percent=round(100 * float(r["blur"]), 4),
               energy_jitter_percent=round(100 * float(r["jit_e"]), 4),
               focusing_jitter_percent=round(100 * float(np.std(r["delta"])), 4),
               equivalent_entry_offset_um=round(float(q[3]), 2),
               screen_dispersion_mid_scan_mm=round(1e3 * float(disp_mid), 4),
               within_shot_skewness=round(float(SKEW[0]), 4),
               pointing_jitter_mrad=round(1e-3 * float(r["sig_th"]), 5),
               spectrometer_pointing_percent=round(100 * float(r["point_rel"]), 4),
               spread_without_pointing_correction_percent=round(100 * float(r["spread_with_pointing"]), 4),
               passive_lens_gradient_T_per_m_at_mean_charge=round(float(r["passive_T_per_m"]), 4),
               reduced_chi2=round(float(r["chi2"]), 4), rounds=int(r["rounds"]))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
