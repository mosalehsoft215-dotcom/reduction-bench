"""Structurally independent analysis of the focus scan (verification evidence).

Written separately from authoring/provenance/lens_model.py and solution/, sharing no code with them.
Like the reference it treats every shot on its own and lets the lens steer the beam; nothing else is
common:
  - the beam is described by its second moments at the lens entrance (emittance, beta, alpha), which
    are fitted and converted to the source at the end, not by the source itself;
  - the momentum distribution within a shot is measured, not assumed Gaussian or of any family: the
    spectrometer shots are each shifted onto the run's mean energy and averaged, and that spectrum is
    deconvolved on the screen (Richardson-Lucy, Gaussian kernel of the beam's own spot on that screen,
    recomputed from this analysis's beam in every round); its bins are the weights;
  - the lens matrix is integrated numerically (RK4 on x'' = -k x, 400 steps, tabulated against k);
  - the steering is a dispersion at the screen, with no ray model: a slice of relative momentum 1 + delta
    lands at -D ln(1 + delta); D is one common number plus a number of every shot's own, for its pointing.
    Near focus the shot's focusing offset and dispersion are scanned together on a grid. Only at the end is
    each shot's dispersion turned into a pointing angle with the lens matrix, and the rms pointing (with every
    shot's fit error) is removed from the spectrometer's run spread through R11 and R12;
  - the data are not shifted: the model is placed on the measured centroid of the shot;
  - the screen is worked out in vectors (where a displaced electron meets the tilted face, how long
    its path through the crystal looks from the camera), the rms of the optics is the Gaussian whose
    edge response has the stated 10-90 % rise, and optics and line of light are one numerical kernel
    on a fine grid;
  - the beam's own wake in the discharge plasma: a shot's charge is its counts in the fit window; the focusing
    offsets are regressed on the set-point trend and on charge / k, and the charge coefficient becomes an extra
    focusing of every momentum bin of the shot in proportion to its charge (none at zero charge); it changes the
    transported sizes only, since the steering is the fitted dispersion;
  - the focusing offset of a shot is found with a coarse scan (wide in the first round, which also reaches the
    wake's own focusing) and a bounded Brent minimiser; the pointing is searched from the second round on.
    python independent_check.py     (task image, run from authoring/evidence)  -> independent_check.txt
Parts: (1) the shape at the true beam against the generator's model; (2) the fit of the task data and
four further draws.
"""
import json
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, least_squares, minimize_scalar
from scipy.special import ndtri

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "environment" / "data"
KEYS = ["source_size_um", "divergence_mrad", "source_distance_mm", "energy_MeV", "energy_spread_percent"]
TOL = {"source_size_um": 3.0, "divergence_mrad": 1.5, "source_distance_mm": 1.5, "energy_MeV": 0.5, "energy_spread_percent": 4.5}
M_E = 0.51099895
LENS_L, DRIFT_L, GRAD_PER_A = 0.2, 1.5, 0.2
CRYSTAL_T, TILT_DEG, FACE_PIXEL, FACE_EDGE = 50e-6, 45.0, 6.80e-6, 43.5e-6    # setup.txt
SPEC_R11, SPEC_R12, SPEC_SCREEN_MM, SPEC_PIXEL_MM = 0.92, 1.85, 0.10, 0.10    # setup.txt
STEP = 2e-6                                                                   # fine grid of the shapes, m


def screen_geometry():
    """the screen in vectors: beam along b, crystal normal n (TILT_DEG from b, crystal turned about the
    vertical), camera along n, so the image shows the face coordinate r . e. Returns the length on the face
    per unit of the beam's horizontal coordinate, and the length on the face of one electron's path"""
    t = np.radians(TILT_DEG)
    b, xh = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    n, e = np.array([np.sin(t), 0.0, np.cos(t)]), np.array([np.cos(t), 0.0, -np.sin(t)])
    hit = xh - b * (n @ xh) / (n @ b)              # where an electron displaced by one unit of x meets the face
    path = b * CRYSTAL_T / (n @ b)                 # from the entrance face to the exit face along the beam
    return abs(hit @ e), abs(path @ e)


FACE_PER_X, STREAK_FACE = screen_geometry()
PIXEL = FACE_PIXEL / FACE_PER_X
RESOLUTION = FACE_EDGE / FACE_PER_X / (ndtri(0.9) - ndtri(0.1))   # Gaussian whose edge rises 10-90 % over FACE_EDGE
DEPTH = STREAK_FACE / FACE_PER_X                                   # line of light of one electron, beam coordinate


def screen_kernel():
    """the optics (Gaussian) convolved with the line of light (uniform), on the fine grid"""
    ng = int(np.ceil(6 * RESOLUTION / STEP))
    g = np.exp(-(np.arange(-ng, ng + 1) * STEP) ** 2 / (2 * RESOLUTION ** 2))
    bx = np.arange(-int(np.ceil(DEPTH / 2 / STEP)) - 1, int(np.ceil(DEPTH / 2 / STEP)) + 2) * STEP
    box = np.clip(np.minimum(bx + STEP / 2, DEPTH / 2) - np.maximum(bx - STEP / 2, -DEPTH / 2), 0, None)
    return np.convolve(g / g.sum(), box / box.sum())


KERNEL = screen_kernel()


def read(path):
    return np.genfromtxt(path, delimiter=",", skip_header=1)


def run_spectrum(d):
    raw = read(d / "spectrometer.csv")[:, 1:]
    energy = read(d / "energy_axis.csv")[:, 1]
    dEdpx = np.gradient(energy)
    px = np.arange(raw.shape[1])
    side = np.r_[px[:50], px[-50:]]
    dens, shot_t = [], []
    for y in raw:
        # Keep signed background-subtracted counts for moment estimates.
        # Clipping every noisy pixel positive biases variance on a wide detector.
        s = (y - np.polyval(np.polyfit(side, y[side], 1), px)) / dEdpx
        s = s / np.sum(s * dEdpx)
        dens.append(s)
        shot_t.append(np.sum(s * dEdpx * energy))
    shot_t = np.array(shot_t)
    w = np.sum(dens, axis=0) * dEdpx
    w /= w.sum()
    t_mean = np.sum(w * energy)
    rel = np.sqrt(np.sum(w * (energy - t_mean) ** 2)) / t_mean          # run spread as measured (with the blur)
    pc = np.sqrt(energy ** 2 + 2 * energy * M_E)
    pc_mean = np.sqrt(t_mean ** 2 + 2 * t_mean * M_E)
    dl = pc / pc_mean - 1.0
    # spectrum within a shot: every shot moved onto the run mean, then averaged
    wc = np.sum([np.interp(energy + (t - t_mean), energy, s, left=0.0, right=0.0) for s, t in zip(dens, shot_t)], axis=0) * dEdpx
    # Only the coadded quadrature density needs nonnegative sampling weights.
    # Total mean and variance above still use signed, unbiased count moments.
    wc = np.maximum(wc, 0)
    wc /= wc.sum()
    i0 = np.argmin(np.abs(energy - t_mean))
    disp = (np.log(energy[i0 + 1]) - np.log(energy[i0 - 1])) / (2 * SPEC_PIXEL_MM)      # d ln T / ds, per mm
    kin = pc_mean ** 2 / ((t_mean + M_E) * t_mean)                                      # dT/T = kin dp/p
    return t_mean, rel, dl, wc, disp, kin


def unblur(w, spot_px, iterations=300):
    """Richardson-Lucy deconvolution of the centred spectrum (weights per pixel) by a Gaussian spot of rms
    spot_px pixels on the screen, inside +-250 pixels of its centroid"""
    px = np.arange(w.size)
    c = np.sum(w * px) / np.sum(w)
    m = np.where(np.abs(px - c) < 250, np.maximum(w, 0.0), 0.0)
    x = np.arange(-int(6 * spot_px) - 1, int(6 * spot_px) + 2)
    ker = np.exp(-x ** 2 / (2 * spot_px ** 2))
    ker /= ker.sum()
    u = m + 1e-15
    for _ in range(iterations):
        u = u * np.convolve(m / np.maximum(np.convolve(u, ker, "same"), 1e-300), ker[::-1], "same")
    return u / u.sum()


def binned(dl, w):
    """merge neighbouring pixels into bins of 0.2 % in momentum (weight-averaged positions)"""
    idx = np.floor(dl / 2e-3).astype(int)
    bins = np.unique(idx)
    wb = np.array([w[idx == b].sum() for b in bins])
    db = np.array([np.sum(w[idx == b] * dl[idx == b]) / max(w[idx == b].sum(), 1e-300) for b in bins])
    keep = wb > 1e-6 * wb.max()
    db, wb = db[keep], wb[keep] / wb[keep].sum()
    return db - np.sum(wb * db), wb


def spectrometer_spot_mm(eps, beta, alpha):
    """rms horizontal spot of the beam on the spectrometer screen (transport R11, R12 from the setup sheet)"""
    gam = (1 + alpha * alpha) / beta
    var = eps * (SPEC_R11 ** 2 * beta - 2 * SPEC_R11 * SPEC_R12 * alpha + SPEC_R12 ** 2 * gam)
    return np.sqrt(var * 1e6 + SPEC_SCREEN_MM ** 2)


def lens_drift(k):
    """numerical 2x2 matrix of the lens (strengths k, any shape) plus drift; returns M11, M12.
    RK4 on (x, x') for the two unit start vectors, all strengths at once."""
    k = np.asarray(k, float)
    n = 400
    h = LENS_L / n
    res = []
    for start in ((1.0, 0.0), (0.0, 1.0)):
        x, xp = np.full(k.shape, start[0]), np.full(k.shape, start[1])
        for _ in range(n):
            a1x, a1p = xp, -k * x
            a2x, a2p = xp + h / 2 * a1p, -k * (x + h / 2 * a1x)
            a3x, a3p = xp + h / 2 * a2p, -k * (x + h / 2 * a2x)
            a4x, a4p = xp + h * a3p, -k * (x + h * a3x)
            x, xp = x + h / 6 * (a1x + 2 * a2x + 2 * a3x + a4x), xp + h / 6 * (a1p + 2 * a2p + 2 * a3p + a4p)
        res.append(x + DRIFT_L * xp)
    return res[0], res[1]


class Table:
    """lens matrix against k, integrated once"""

    def __init__(self, k_lo, k_hi):
        self.k = np.linspace(k_lo, k_hi, 4001)
        self.m11, self.m12 = lens_drift(self.k)

    def __call__(self, k):
        return np.interp(k, self.k, self.m11), np.interp(k, self.k, self.m12)


def shape(x, m11, m12, cen, eps, beta, alpha, weights):
    """profile at positions x (from the shot centroid): momentum bins with their transported sizes and
    landing places, summed on a fine grid, convolved with the resolution, read at x"""
    s2 = eps * (m11 ** 2 * beta - 2 * m11 * m12 * alpha + m12 ** 2 * (1 + alpha ** 2) / beta)
    c = cen - np.sum(weights * cen)
    half = np.max(np.abs(x)) + 6 * RESOLUTION + DEPTH
    fine = np.arange(-half, half + STEP, STEP)
    f = np.sum(weights[:, None] * np.exp(-(fine[None, :] - c[:, None]) ** 2 / (2 * s2[:, None]))
               / np.sqrt(2 * np.pi * s2[:, None]), axis=0)
    return np.interp(x, fine, np.convolve(f, KERNEL, mode="same"))


def scan(d):
    """every shot: set current, positions from its centroid, counts above background, weights"""
    a = read(d / "focus_scan.csv")
    cur, y = a[:, 0], a[:, 2:]
    n = y.shape[1]
    px = np.arange(n)
    xc = (px - n / 2 + 0.5) * PIXEL
    side = np.r_[px[:150], px[-150:]]
    out = []
    for I, yy in zip(cur, y):
        z = yy - np.polyval(np.polyfit(side, yy[side], 1), px)
        read_noise = np.std(z[side])
        g = lambda t, a_, m, w: a_ * np.exp(-(t - m) ** 2 / (2 * w ** 2))
        q, _ = curve_fit(g, xc, z, p0=[z.max(), xc[np.argmax(z)], 60e-6])
        win = np.abs(xc - q[1]) < np.clip(7 * abs(q[2]), 60 * PIXEL, 420 * PIXEL)
        mid = np.sum(z[win] * xc[win]) / np.sum(z[win])                  # centroid of the shot
        sm = np.convolve(z, np.ones(7) / 7, "same")
        w = 1.0 / np.sqrt(read_noise ** 2 + np.clip(sm[win], 0, None))
        out.append((I, xc[win] - mid, z[win], w, float(np.sum(z[win]))))
    return out


def fit_shot(x, z, w, sh):
    """amplitude and offset of the model shape by weighted least squares; returns the weighted residuals"""
    A = np.vstack([sh * w, w]).T
    return z * w - A @ np.linalg.lstsq(A, z * w, rcond=None)[0]


def analyse(d):
    d = Path(d)
    t_mean, rel_meas, dl_meas, w_shot, disp, kin = run_spectrum(d)
    pc = np.sqrt(t_mean ** 2 + 2 * t_mean * M_E) / 1000.0
    data = scan(d)
    cur = np.array([s[0] for s in data])
    # the discharge plasma: the beam's own wake focuses it in proportion to its charge (counts of the shot),
    # about its own centroid, so it changes the transported sizes and not the screen dispersion
    charge = np.array([s[4] for s in data])
    qrel = charge / charge.mean()
    kwake = [0.0]
    k_nom = 0.299792458 * GRAD_PER_A * cur / pc
    table = Table(0.85 * k_nom.min(), 1.2 * k_nom.max())
    rms_shot = np.sqrt(np.sum(w_shot * dl_meas ** 2) - np.sum(w_shot * dl_meas) ** 2)
    trend = np.vstack([np.ones(cur.size), np.log(cur) - np.log(cur).mean()]).T
    # an error e in the wake's strength leaves the shot's momentum offset at about -e qrel / k_nom
    trend_q = np.vstack([trend.T, qrel / k_nom]).T
    no_wake = bool(__import__('os').environ.get('INDEP_NO_WAKE'))

    def remove_trend(offs):
        if no_wake:
            return offs - trend @ np.linalg.lstsq(trend, offs, rcond=None)[0]
        c = np.linalg.lstsq(trend_q, offs, rcond=None)[0]
        kwake[0] -= c[2]
        return offs - trend_q @ c
    q = np.array([3.0, 2.0, -40.0, 0.5])                     # emittance (nm), beta (m), alpha, dispersion (mm)
    offs = np.zeros(cur.size)
    shrink = 1.0

    extra = np.zeros(cur.size)            # each shot's own screen dispersion from its pointing (mm), mean in q[3]
    extra_se = np.full(cur.size, np.inf)

    def shot_res(i, off, q, deltas, weights, dd=None):
        I, x, z, w, _ = data[i]
        mom = 1 + off + deltas
        m11, m12 = table(k_nom[i] / mom + kwake[0] * qrel[i] / (1 + deltas))
        disp_i = q[3] + (extra[i] if dd is None else dd)
        return fit_shot(x, z, w, shape(x, m11, m12, -disp_i * 1e-3 * np.log(mom), q[0] * 1e-9, q[1], q[2], weights))

    def refine(i, q, deltas, weights, o0, d0):
        g = least_squares(lambda v: shot_res(i, v[0], q, deltas, weights, v[1]), [o0, d0], x_scale=[1e-3, 0.2],
                          diff_step=[1e-6, 1e-4], max_nfev=60)
        try:
            se = float(np.sqrt(np.linalg.inv(g.jac.T @ g.jac)[1, 1]))
        except np.linalg.LinAlgError:
            se = np.inf
        return g.x[0], g.x[1], se

    spot_px = spectrometer_spot_mm(q[0] * 1e-9, q[1], q[2]) / SPEC_PIXEL_MM
    for rnd in range(13):
        # momentum bins of a shot: the centred spectrum, deconvolved on the screen by the beam's own spot
        deltas, weights = binned(dl_meas, unblur(w_shot, spot_px))
        if rnd > 0:                                                       # round 0: all shots at their set current
            for i in range(cur.size):
                cost = lambda o: np.sum(shot_res(i, o, q, deltas, weights) ** 2)
                # round 1 reaches past the wake's own focusing (a few % of k at a heavy shot); later rounds stay close
                grid = offs[i] + (np.linspace(-0.06, 0.06, 31) if rnd == 1 else np.linspace(-0.002, 0.002, 5))
                o0 = grid[int(np.argmin([cost(o) for o in grid]))]
                h = grid[1] - grid[0]
                offs[i] = minimize_scalar(cost, bounds=(o0 - h, o0 + h), method="bounded", options={"xatol": 2e-5}).x
            offs = remove_trend(offs)   # random jitter: no mean, no trend; the charge part is the wake's
            if rnd > 1:                 # the pointing is searched once the wake has a first value
                # pointing: a dispersion of its own for every shot. Near focus a shot is also matched with its focus
                # on the other side and the dispersion reversed, so offset and dispersion are scanned together there.
                for i in range(cur.size):
                    cost = lambda o, dd: np.sum(shot_res(i, o, q, deltas, weights, dd) ** 2)
                    if data[i][1].size <= 400:
                        cand = [(o, dd) for o in offs[i] + np.linspace(-0.02, 0.02, 21) for dd in np.linspace(-4.5, 4.5, 25)]
                    elif rnd <= 3:
                        cand = [(offs[i], dd) for dd in np.linspace(-4.5, 4.5, 25)]
                    else:
                        cand = [(offs[i], extra[i])]
                    o0, d0 = cand[int(np.argmin([cost(o, dd) for o, dd in cand]))]
                    offs[i], extra[i], extra_se[i] = refine(i, q, deltas, weights, o0, d0)
                offs = remove_trend(offs)
                extra -= extra.mean()
        f = least_squares(lambda v: np.concatenate([shot_res(i, offs[i], v, deltas, weights) for i in range(cur.size)]),
                          q, bounds=([1e-3, 0.01, -500.0, -20.0], [1e3, 100.0, 500.0, 20.0]), x_scale="jac", max_nfev=40)
        change = np.max(np.abs(f.x[:3] / q[:3] - 1))
        q = f.x
        spot_px = spectrometer_spot_mm(q[0] * 1e-9, q[1], q[2]) / SPEC_PIXEL_MM
        blur = spectrometer_spot_mm(q[0] * 1e-9, q[1], q[2]) * disp
        shrink = np.sqrt(max(1.0 - (blur / rel_meas) ** 2, 1e-6))
        if change < 3e-4 and rnd > 3:
            break
    eps, beta, alpha = q[0] * 1e-9, q[1], q[2]
    gam_tw = (1 + alpha ** 2) / beta                                     # second moments -> waist of the beam
    size_um, div_mrad, dist_mm = 1e6 * np.sqrt(eps / gam_tw), 1e3 * np.sqrt(eps * gam_tw), -1e3 * alpha / gam_tw
    dof = f.fun.size - 4 - 4 * cur.size
    chi2 = 2 * f.cost / dof
    # dispersion per unit pointing: a ray leaving the source at angle th crosses the lens entrance at dist*th
    # with angle th, and lands at (M11 dist + M12) th; the slice of momentum 1+delta has k/(1+delta)
    keff = k_nom / (1 + offs)
    dk = 1e-4 * keff
    lever = lambda k: table(k)[0] * dist_mm * 1e-3 + table(k)[1]
    per_rad = keff * (lever(keff + dk) - lever(keff - dk)) / (2 * dk) * 1e3                # mm of dispersion per rad
    th, th_se = extra / per_rad, extra_se * np.sqrt(max(chi2, 1.0)) / np.abs(per_rad)
    th = th - th.mean()
    nll = lambda s: np.sum(np.log(s * s + th_se ** 2) + th ** 2 / (s * s + th_se ** 2))
    sig_th = minimize_scalar(nll, bounds=(1e-6, 5e-3), method="bounded").x
    # the spectrometer moves a shot by (R11 dist + R12) th: remove that from the run spread
    point_rel = (SPEC_R11 * dist_mm * 1e-3 + SPEC_R12) * sig_th * 1e3 * disp
    spread = 100 * np.sqrt(max((rel_meas * shrink) ** 2 - point_rel ** 2, 1e-12))
    return [size_um, div_mrad, dist_mm, t_mean, spread], chi2, q[3], np.std(offs), 1e3 * sig_th, kwake[0] * pc / 0.299792458


def _job(seed):
    if seed is None:
        return "task data", analyse(DATA)
    with tempfile.TemporaryDirectory() as d:
        subprocess.run([sys.executable, str(HERE.parent / "provenance" / "generate_focus_scan.py"), str(seed), d], check=True,
                       cwd=HERE.parent / "provenance")
        return f"noise draw {seed}", analyse(d)


def forward_check():
    sys.path.insert(0, str(HERE.parent / "provenance"))
    import generate_focus_scan as G
    from lens_model import profile
    pcm = np.sqrt(G.T_MEAN ** 2 + 2 * G.T_MEAN * M_E) / 1000.0
    from lens_model import skew_nodes
    xg, wg = skew_nodes(G.SKEW_SHOT)                 # the generator's momentum distribution within a shot
    xc = (np.arange(1024) - 512 + 0.5) * PIXEL
    worst = 0.0
    for I in G.scan_currents()[::5]:
        k = 0.299792458 * GRAD_PER_A * I / pcm
        m11, m12 = lens_drift(k / (1 + xg * G.SPREAD_SHOT))
        sel = np.abs(xc) < 2.0e-3
        mine = shape(xc[sel], m11, m12, m11 * G.ENTRY[0] + m12 * G.ENTRY[1], G.EPS, G.BETA, G.ALPHA, wg)
        theirs = profile(xc[sel], k, G.SIZE, G.DIV, G.DIST, G.SPREAD_SHOT, G.RES, entry=G.ENTRY, box=G.BOX, nodes=(xg, wg))
        worst = max(worst, np.max(np.abs(mine - theirs)) / theirs.max())
    return [f"screen from its geometry in vectors: {PIXEL * 1e6:.4f} um per pixel, optics {RESOLUTION * 1e6:.4f} um rms, line of "
            f"light {DEPTH * 1e6:.4f} um in the beam's coordinate (generator: {G.PIX * 1e6:.4f}, {G.RES * 1e6:.4f}, {G.BOX * 1e6:.4f})",
            f"profile shape of a shot at the true beam, every 5th current: max |difference| {100 * worst:.3f} % of the peak"]


if __name__ == "__main__":
    truth = json.loads((ROOT / "tests" / "ground_truth.json").read_text())
    out = ["Structurally independent analysis (second moments at the lens entrance, deconvolved measured spectrum as weights, RK4 lens "
           "matrix, steering as a screen dispersion, model placed on the shot centroid, screen geometry in vectors with a numerical kernel)", "",
           "(1) forward check against the generator:"]
    out += forward_check()
    out += ["", "(2) fit with the independent analysis (every shot with its own focusing offset):"]
    t0 = time.time()
    with Pool(5) as pool:
        res = pool.map(_job, [None, 301, 302, 303, 304])
    vals = dict(res)["task data"][0]
    task = dict(values=dict(zip(KEYS, map(float, vals))), errors_percent={k: 100 * (float(v) / truth[k] - 1) for k, v in zip(KEYS, vals)},
                reduced_chi2=float(dict(res)["task data"][1]), pointing_mrad=float(dict(res)["task data"][4]),
                runtime_sec=time.time() - t0, note="Independent second-moment/RK4 source reduction with a Richardson-Lucy "
                "deconvolved spectrum and a dispersion of every shot's own for its pointing; task data (run alongside four noise draws).")
    (HERE / "independent_source_repaired.json").write_text(json.dumps(task, indent=2) + chr(10), encoding="utf-8", newline=chr(10))
    worst = {k: 0.0 for k in KEYS}
    for name, (vals, chi2, disp_mm, jit, sig_mrad, g_wake) in res:
        err = {k: 100 * (vals[i] / truth[k] - 1) for i, k in enumerate(KEYS)}
        for k in KEYS:
            worst[k] = max(worst[k], abs(err[k]))
        out.append(f"{name:14s} " + "  ".join(f"{k} {err[k]:+.2f}%" for k in KEYS)
                   + f"  reduced chi2 {chi2:.3f}  dispersion {disp_mm:+.3f} mm  focusing jitter {100 * jit:.3f} %  pointing {sig_mrad:.3f} mrad"
                   + f"  wake {g_wake:.3f} T/m at mean charge")
    out.append("")
    out.append("largest |error| of the independent analysis: " + ", ".join(
        f"{k} {worst[k]:.2f} % (bound {TOL[k]} %, {TOL[k] / max(worst[k], 1e-9):.1f}x inside)" for k in KEYS))
    (HERE / "independent_check.txt").write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print("\n".join(out))
