"""An independently written correct analysis. It shares no code with solution/ and makes a
different, equally legitimate choice at every step:

  step                 reference (solution/wlfit.py)           this script
  NFW lensing          Wright & Brainerd (2000), wlsim.py       Bartelmann (1996) closed forms, own cosmology
  field template       calibration z_phot, Gaussian smoothing   calibration z_phot, triangular kernel
  member photo-z model unblended core; blended: shifted part     the same components plus the upper end z_max of
                       (half-normal) + uniform outliers to 2.5  the outliers as a free parameter
  shape fit            groups of 15 clusters, inner 6 bins,     6 (3, 12) redshift slices, inner 6 bins, field,
                       3 counts per group-bin (Newton)          core and blended counts per slice-bin (Poisson)
  member fractions     per group-bin                            per slice and bin
  shear fit            6000 cluster-bin points, per-point       shear stacked per slice and bin (weights
                       dilution                                 n_src), dilution of the stacked model

It also repeats the fit with misspecified member models (one blended share for all bins, which
is what an analysis that ignores the crowding assumes; no shifted part; outliers spread to
z = 3; a symmetric instead of a one-sided shifted part) and with a different one-sided shape
(exponential), reports how strongly the data prefer the fitted model, and fits the member
shapes on the inner 4 or all 10 bins instead of the inner 6.
usage: python independent_correct.py DATA_DIR
"""
import math, os, sys
import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import ndtr

TRUTH = dict(M=7.60, c=6.80)
TOL = dict(M=0.03, c=0.09)
CUT = 0.1                                       # source cut: z_phot > z_l + CUT
C_KMS, G, H0, OM = 299792.458, 4.30091727e-9, 70.0, 0.3

_zg = np.linspace(0.0, 5.0, 20001)
_ie = 1 / np.sqrt(OM * (1 + _zg) ** 3 + 1 - OM)
_dc = np.r_[0.0, np.cumsum(0.5 * (_ie[1:] + _ie[:-1]) * np.diff(_zg))] * C_KMS / H0


def comoving(z):
    return np.interp(z, _zg, _dc)


def sigma_crit_inv(zl, zs):
    """1 / Sigma_crit in Mpc^2 / Msun (flat universe); zero for sources in front of the lens."""
    with np.errstate(divide="ignore", invalid="ignore"):
        x = 4 * math.pi * G / C_KMS ** 2 * comoving(zl) / (1 + zl) * (1 - comoving(zl) / comoving(zs))
    return np.where(zs > zl, x, 0.0)


def rho_crit(z):
    return 3 * H0 ** 2 * (OM * (1 + z) ** 3 + 1 - OM) / (8 * math.pi * G)


def nfw_lensing(R, M, c, zl):
    """Sigma and Delta Sigma (Msun / Mpc^2) of an NFW halo (Bartelmann 1996).
    R: (n_clusters, n_radii) physical Mpc; zl: (n_clusters,)."""
    rc = rho_crit(zl)[:, None]
    rs = (3 * M / (800 * math.pi * rc)) ** (1 / 3) / c
    ks = 2 * rs * 200 / 3 * c ** 3 / (math.log(1 + c) - c / (1 + c)) * rc      # 2 rho_s r_s
    x = R / rs
    x = np.where(np.abs(x - 1) < 1e-6, 1 + 1e-6, x)
    F = np.where(x < 1, np.arccosh(1 / np.minimum(x, 1)) / np.sqrt(np.abs(1 - x * x)),
                 np.arccos(1 / np.maximum(x, 1)) / np.sqrt(np.abs(x * x - 1)))
    sig = ks * (1 - F) / (x * x - 1)
    mean_inside = 2 * ks * (F + np.log(x / 2)) / (x * x)
    return sig, mean_inside - sig


def load(d):
    rd = lambda f: np.genfromtxt(os.path.join(d, f), delimiter=",", names=True)
    cl, P, cal = rd("clusters.csv"), rd("shear_profiles.csv"), rd("photoz_calibration.csv")
    with open(os.path.join(d, "zphot_histograms.csv")) as fh:
        head = fh.readline().strip().split(",")
    H = np.loadtxt(os.path.join(d, "zphot_histograms.csv"), delimiter=",", skiprows=1)
    zlo = np.array([float(h[2:]) for h in head[2:]])
    ncl, nb = len(cl), int(P["bin"].max()) + 1
    zl = np.empty(ncl); zl[cl["cluster"].astype(int)] = cl["z"]
    k, b = P["cluster"].astype(int), P["bin"].astype(int)

    def grid(col):
        a = np.zeros((ncl, nb)); a[k, b] = np.nan_to_num(P[col]); return a
    hist = np.zeros((ncl, nb, len(zlo)))
    hist[H[:, 0].astype(int), H[:, 1].astype(int)] = H[:, 2:]
    return dict(zl=zl, n=grid("n_src"), et=grid("et_mean"), err=grid("e_err"), hist=hist,
                zedges=np.r_[zlo, 2 * zlo[-1] - zlo[-2]], edges=np.unique(np.r_[P["R_lo_Mpc"], P["R_hi_Mpc"]]),
                zs=cal["z_spec_lo"] + 0.005, zp=cal["z_phot_lo"], ncal=cal["n"])


def field_template(D):
    """Field z_phot distribution per histogram cell: the calibration counts summed over true
    redshift (the calibration cells in z_phot are the histogram cells)."""
    t = np.histogram(D["zp"] + 0.01, D["zedges"], weights=D["ncal"])[0]
    t = np.maximum(t, 1e-3)
    return t / t.sum()


def shifted_cdf(x, s2, kind):
    """CDF of the shifted member component in x = (z_phot - z_l) / (1 + z_l): one-sided
    half-normal (the default), one-sided exponential, or a symmetric Gaussian."""
    if kind == "exp":
        return np.where(x > 0, 1 - np.exp(-np.maximum(x, 0) / s2), 0.0)
    if kind == "gauss":
        return ndtr(x / s2)
    return np.clip(2 * ndtr(x / s2) - 1, 0, None)


def core_cells(zl, zedges, s1):
    """Photo-z core of the unblended members per histogram cell: Gaussian of width s1 (1 + z_l)
    around z_l, photo-z below 0 recorded as 0; one row per cluster redshift."""
    cdf = ndtr((zedges[None, :] - zl[:, None]) / (s1 * (1 + zl[:, None])))
    core = np.diff(cdf, axis=1)
    core[:, 0] += cdf[:, 0]
    return core / core.sum(1, keepdims=True)


def blend_cells(zl, zedges, s2, share, zmax, kind="half"):
    """Photo-z of the blended members per histogram cell: a share pulled to higher photo-z
    (shifted_cdf, width s2) and the rest catastrophic, uniform in [0, zmax]."""
    cs = shifted_cdf((zedges[None, :] - zl[:, None]) / (1 + zl[:, None]), s2, kind)
    shifted = np.diff(cs, axis=1)
    shifted[:, 0] += cs[:, 0]
    out = np.clip(np.minimum(zedges[1:], zmax) - zedges[:-1], 0, None) / zmax
    q = share * shifted + (1 - share) * out[None, :]
    return q / q.sum(1, keepdims=True)


def sel_core(zl, s1):
    """Share of the unblended members above the cut z_phot > z_l + CUT."""
    return ndtr(-CUT / (s1 * (1 + zl)))


def sel_blend(zl, zedges, s2, share, zmax, kind="half"):
    """Selected blended members (beyond the histogram range included) per blended member
    counted in the histograms."""
    top, z1 = zedges[-1], 1 + zl
    sel = share * (1 - shifted_cdf(CUT / z1, s2, kind)) + (1 - share) * np.clip(zmax - (zl + CUT), 0, None) / zmax
    inside = share * shifted_cdf((top - zl) / z1, s2, kind) + (1 - share) * min(top, zmax) / zmax
    return sel / inside


def amplitudes(h, T, pF, iters=40):
    """Poisson maximum-likelihood field and member counts of every row of h, for the member
    templates T (one per row) and the field template pF: EM steps, then Newton steps on the
    concave log-likelihood kept inside A, B > 0."""
    P = np.broadcast_to(pF, h.shape)
    A = 0.5 * h.sum(1) + 1.0; B = A.copy()
    for it in range(iters):
        mu = A[:, None] * P + B[:, None] * T
        if it < 10:
            A, B = A * (h * P / mu).sum(1), B * (h * T / mu).sum(1)
            continue
        r = h / mu; r2 = r / mu
        gA, gB = (r * P).sum(1) - 1, (r * T).sum(1) - 1
        a, ab, bb = (r2 * P * P).sum(1), (r2 * P * T).sum(1), (r2 * T * T).sum(1)
        det = a * bb - ab * ab
        dA, dB = (bb * gA - ab * gB) / det, (a * gB - ab * gA) / det
        t = np.ones_like(A)
        for X, dX in ((A, dA), (B, dB)):
            m = X + dX <= 0
            t[m] = np.minimum(t[m], 0.5 * X[m] / -dX[m])
        A, B = A + t * dA, B + t * dB
    mu = A[:, None] * P + B[:, None] * T
    return A, B, float((h * np.log(mu + 1e-300) - mu).sum())


def amplitudes3(h, C, Bl, pF, iters=60):
    """Poisson maximum-likelihood counts of field, unblended and blended members in every row
    of h (templates pF, C and Bl): multiplicative EM steps, then Newton steps."""
    T = np.stack([np.broadcast_to(pF, h.shape), C, Bl], axis=1)
    A = np.repeat(h.sum(1, keepdims=True) / 3.0, 3, axis=1) + 1.0
    for it in range(iters):
        mu = np.einsum("nk,nkc->nc", A, T) + 1e-300
        if it < 15:
            A = A * np.einsum("nc,nkc->nk", h / mu, T)
            continue
        r = h / mu
        g = np.einsum("nc,nkc->nk", r, T) - 1.0
        Hn = np.einsum("nc,nkc,njc->nkj", r / mu, T, T) + 1e-9 * np.eye(3)
        st = np.linalg.solve(Hn, g[..., None])[..., 0]
        lim = np.where(A + st <= 0, 0.5 * A / np.maximum(-st, 1e-300), np.inf).min(1)
        A = A + np.minimum(1.0, lim)[:, None] * st
    mu = np.einsum("nk,nkc->nc", A, T) + 1e-300
    return A, float((h * np.log(mu) - mu).sum())


def slices(zl, n):
    return [np.sort(g) for g in np.array_split(np.argsort(zl), n)]


PARAMS = ("s1", "s2", "share", "zmax")
BOUNDS = dict(s1=(0.01, 0.4), s2=(0.02, 1.5), share=(0.0, 1.0), zmax=(0.6, 8.0), P=(0.0, 0.99))


def _templates(D, groups, d, kind):
    C = np.array([core_cells(D["zl"][g], D["zedges"], d["s1"]).mean(0) for g in groups])
    Bl = np.array([blend_cells(D["zl"][g], D["zedges"], d["s2"], d["share"], d["zmax"], kind).mean(0) for g in groups])
    return C, Bl


def fit_member_model(D, groups, pF, fixed=None, bins=range(6), kind="half", constant=False):
    """Member photo-z shapes fitted to the histograms of the inner bins stacked per redshift
    slice and radial bin. Blending follows the crowding, so every slice-bin has its own field,
    unblended-member and blended-member counts. constant=True instead imposes one blended share
    P on all slice-bins (one member photo-z mixture everywhere). Returns the parameters and
    -ln L."""
    fixed = dict(fixed or {})
    names = PARAMS + (("P",) if constant else ())
    free = [k for k in names if k not in fixed]
    nb = len(bins)
    h = np.array([D["hist"][g, b].sum(0) for g in groups for b in bins])

    def unpack(x):
        d = dict(fixed); d.update(zip(free, x)); return d

    def nll(x):
        d = unpack(x)
        if any(not BOUNDS[k][0] <= d[k] <= BOUNDS[k][1] for k in free):
            return 1e30
        C, Bl = _templates(D, groups, d, kind)
        C, Bl = np.repeat(C, nb, axis=0), np.repeat(Bl, nb, axis=0)
        if constant:
            return -amplitudes(h, (1 - d["P"]) * C + d["P"] * Bl, pF)[2]
        return -amplitudes3(h, C, Bl, pF)[1]
    best = None
    for sh0 in (0.4, 0.8):
        for zm0 in (2.0, 3.0):
            st = dict(s1=0.06, s2=0.3, share=sh0, zmax=zm0, P=0.4)
            r = minimize(nll, [st[k] for k in free], method="Nelder-Mead",
                         options=dict(maxiter=3000, maxfev=3000, xatol=1e-5, fatol=1e-4))
            if best is None or r.fun < best.fun:
                best = r
    out = unpack(best.x)
    out["kind"] = kind
    return out, float(best.fun)


def member_fractions(D, groups, pF, shape):
    """Member fraction of the selected sources per slice and radial bin, and the blended share
    of the members per radial bin."""
    nb = D["hist"].shape[1]
    f = np.zeros((len(groups), nb)); Ac = np.zeros(nb); Ab = np.zeros(nb)
    kw = dict(s2=shape["s2"], share=shape["share"], zmax=shape["zmax"], kind=shape["kind"])
    for i, g in enumerate(groups):
        zl = D["zl"][g]
        C = core_cells(zl, D["zedges"], shape["s1"]).mean(0)
        Bl = blend_cells(zl, D["zedges"], **kw).mean(0)
        h = D["hist"][g].sum(0)
        pc, pb = sel_core(zl, shape["s1"]).mean(), sel_blend(zl, D["zedges"], **kw).mean()
        if "P" in shape:
            _, B, _ = amplitudes(h, np.repeat(((1 - shape["P"]) * C + shape["P"] * Bl)[None, :], nb, 0), pF)
            c_, b_ = (1 - shape["P"]) * B, shape["P"] * B
        else:
            A, _ = amplitudes3(h, np.repeat(C[None, :], nb, 0), np.repeat(Bl[None, :], nb, 0), pF)
            c_, b_ = A[:, 1], A[:, 2]
        f[i] = (c_ * pc + b_ * pb) / D["n"][g].sum(0)
        Ac += c_; Ab += b_
    return f, Ab / (Ac + Ab)


class Lens:
    """Reduced-shear NFW model per cluster and radial bin, averaged over the true redshifts of
    the calibration galaxies that pass each cluster's cut and over the bin area."""
    def __init__(self, D, nodes=4, dz=0.025):
        self.zl = D["zl"]
        ze = np.arange(0.0, D["zs"].max() + dz, dz)
        dzp = D["zedges"][1] - D["zedges"][0]
        # calibration galaxies above each cluster's cut: the share of a z_phot cell above it
        w = np.array([np.histogram(D["zs"], ze, weights=D["ncal"] * np.clip((D["zp"] + dzp - z - CUT) / dzp, 0, 1))[0]
                      for z in self.zl])
        used = w.sum(0) > 0
        self.w = w[:, used] / w.sum(1, keepdims=True)
        self.S = sigma_crit_inv(self.zl[:, None], (0.5 * (ze[1:] + ze[:-1]))[None, used])
        x, wt = np.polynomial.legendre.leggauss(nodes)
        u0, u1 = D["edges"][:-1] ** 2, D["edges"][1:] ** 2                     # nodes uniform in area
        self.R = np.sqrt(0.5 * (u1 - u0)[:, None] * x + 0.5 * (u1 + u0)[:, None]).ravel()
        self.wt, self.nodes = wt / wt.sum(), nodes

    def bins(self, M, c):
        sig, ds = nfw_lensing(np.broadcast_to(self.R, (len(self.zl), self.R.size)), M, c, self.zl)
        S = self.S[:, :, None]
        g = np.einsum("kz,kzr->kr", self.w, S * ds[:, None, :] / (1 - S * sig[:, None, :]))
        return (g.reshape(len(self.zl), -1, self.nodes) * self.wt).sum(2)


def fit_shear(D, groups, f, lens, rmin=0.0):
    use = D["edges"][:-1] >= rmin - 1e-9
    W = [D["n"][g] for g in groups]                                             # weights n_src
    ebar = np.array([(w * D["et"][g]).sum(0) / w.sum(0) for w, g in zip(W, groups)])
    sbar = np.array([np.sqrt((w ** 2 * D["err"][g] ** 2).sum(0)) / w.sum(0) for w, g in zip(W, groups)])

    def resid(p):
        gb = lens.bins(10 ** p[0], p[1])
        mod = np.array([(1 - f[i]) * (w * gb[g]).sum(0) / w.sum(0) for i, (w, g) in enumerate(zip(W, groups))])
        return ((mod - ebar) / sbar)[:, use].ravel()
    r = least_squares(resid, [14.8, 4.0], bounds=([13.5, 0.5], [16.0, 20.0]), x_scale=[0.05, 0.5])
    return 10 ** r.x[0] / 1e14, float(r.x[1])


def report(label, M, c):
    ok = abs(M / TRUTH["M"] - 1) <= TOL["M"] and abs(c / TRUTH["c"] - 1) <= TOL["c"]
    print("%-62s M200c %.3f (%+5.1f %%)  c200c %.3f (%+5.1f %%)  %s" % (
        label, M, 100 * (M / TRUTH["M"] - 1), c, 100 * (c / TRUTH["c"] - 1), "pass" if ok else "fail"), flush=True)


def describe(sh):
    return ("core %.4f, shifted width %.3f, shifted share of the blended %.3f, outliers up to z = %.2f"
            % (sh["s1"], sh["s2"], sh["share"], sh["zmax"])
            + (", one blended share %.3f for all bins" % sh["P"] if "P" in sh else "")
            + ("" if sh.get("kind", "half") == "half" else " (%s)" % sh["kind"]))


def run(d):
    D = load(d)
    pF = field_template(D)
    lens = Lens(D)
    for nsl in (6, 3, 12):
        groups = slices(D["zl"], nsl)
        shape, nll = fit_member_model(D, groups, pF)
        f, bl = member_fractions(D, groups, pF, shape)
        print("%d slices: %s" % (nsl, describe(shape)))
        if nsl == 6:
            N = np.array([D["n"][g].sum(0) for g in groups])
            print("blended share of the members: " + " ".join("%.3f" % x for x in bl))
            print("member fraction of the sources: " + " ".join("%.3f" % x for x in (f * N).sum(0) / N.sum(0)))
        report("independent, %d redshift slices, all bins" % nsl, *fit_shear(D, groups, f, lens))
        if nsl != 6:
            continue
        report("independent, 6 redshift slices, R > 0.5 Mpc", *fit_shear(D, groups, f, lens, rmin=0.5))
        for nbins in (4, 10):
            sh2, _ = fit_member_model(D, groups, pF, bins=range(nbins))
            report("independent, member shapes fitted on bins 0-%d" % (nbins - 1),
                   *fit_shear(D, groups, member_fractions(D, groups, pF, sh2)[0], lens))
        for lab, kw in (("misspecified: one member photo-z mixture in all bins", dict(constant=True)),
                        ("misspecified: no shifted part (blended members uniform only)", dict(fixed=dict(share=0.0, s2=0.3))),
                        ("misspecified: outliers spread to z = 3", dict(fixed=dict(zmax=3.0))),
                        ("independent, shifted part one-sided exponential", dict(kind="exp")),
                        ("misspecified: shifted part symmetric (broad Gaussian)", dict(kind="gauss"))):
            sh2, nll2 = fit_member_model(D, groups, pF, **kw)
            print("%s: %s; -2 Delta ln L = %.0f" % (lab.split(": ")[-1].replace("independent, ", ""), describe(sh2), 2 * (nll2 - nll)))
            report(lab, *fit_shear(D, groups, member_fractions(D, groups, pF, sh2)[0], lens))


if __name__ == "__main__":
    for d in sys.argv[1:]:
        print("==", d, flush=True)
        run(d)
