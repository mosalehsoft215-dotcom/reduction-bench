"""Reference analysis: phase-connected timing solution of a binary pulsar from sparse, noisy sessions.

    python psr_fit.py <data dir> <results.json>

1. Turn counts inside a session. A search over apparent frequency and drift gives several ways of counting the
   turns between the arrival times of one session that fit about equally well; all of them are kept.
2. Orbit. The apparent frequency and drift of every session must lie on the Doppler curve of one elliptical orbit. For each
   trial orbital frequency a harmonic series is fitted, every session is moved to the counting that lies nearest
   to the curve, and the fit is repeated; the best orbital periods go on to a fit of the full orbit.
3. Turn counts between sessions. Around the orbit of step 2 the phase is linear in the parameters, so the
   whole numbers of turns between sessions solve an integer least-squares problem: the lattice is reduced
   (Lenstra-Lenstra-Lovasz), every candidate close to the best is listed, and each is refitted with the full
   model. The one with the smallest chi-square over all arrival times is the solution.
"""
import json
import sys
from pathlib import Path

import numpy as np

from psr_model import ARCSEC, DAY, E_EAST, E_NORTH, T_REF, emission, phase, roemer_rate

F_APPARENT = 287.46                                                  # Hz, timing_notes.txt

STEPS = np.array([0, 0, 0, 0, 0, 1e-2, 1e-5, 1e-6, 1e-5, 1e-1])     # finite-difference steps of the orbit
QSTEP = np.array([1e-7, 1e-1, 1e-5, 1e-5, 1e-4, 1.0, 1e-2, 1e-2])


def sessions(t):
    cut = np.r_[0, np.where(np.diff(t) > 0.2 * DAY)[0] + 1, len(t)]
    return [np.arange(a, b) for a, b in zip(cut[:-1], cut[1:])]


def local(t, f_mid, half=0.12, drift=6e-6, keep=12, peaks=60):
    """ways of counting the turns inside one session: list of (rms, counts, f, fdot), best first"""
    tau = t - t.mean()
    span = np.ptp(tau)
    f = np.arange(-half, half, 0.08 / span)
    fd = np.arange(-drift, drift, 0.4 / span ** 2)
    z = np.array([np.abs(np.exp(2j * np.pi * (np.outer(f_mid + f, tau) + 0.5 * d * tau ** 2)).sum(1)) for d in fd])
    a = np.c_[np.ones_like(tau), tau, tau ** 2]
    found = {}
    for flat in np.argsort(z, axis=None)[::-1][:peaks * 40]:
        i, j = np.unravel_index(flat, z.shape)
        c = np.array([0.0, f_mid + f[j], fd[i] / 2])
        for _ in range(4):
            pred = a @ c
            n = np.round(pred - pred[0])
            c = np.linalg.lstsq(a, n, rcond=None)[0]
        key = tuple(n.astype(int))
        if key not in found:
            found[key] = (float(np.sqrt(np.mean((a @ c - n) ** 2))), n, c[1], 2 * c[2])
        if len(found) >= peaks:
            break
    out = sorted(found.values(), key=lambda q: q[0])
    return [q for q in out if q[0] < 2.5 * out[0][0] + 0.01 and abs(q[3]) < 1.5 * drift][:keep]


def harmonic_scan(ts, cands, rdot, harmonics=3, weight=600.0, rounds=3, keep=4, most=6):
    """orbital periods at which one counting per session lies on a smooth periodic Doppler curve.

    For every trial orbital frequency a harmonic series (plus the Doppler shift of a position offset) is fitted to
    the apparent frequencies and drifts; the fit is repeated on the best four fifths of the sessions and every
    session is moved to the counting nearest the curve. Halves and thirds of the best periods are added, because a
    series at twice the true period fits as well."""
    n_s, width = len(ts), max(len(c) for c in cands)
    freq = np.full((n_s, width), np.inf)
    drift = np.full((n_s, width), np.inf)
    for s, c in enumerate(cands):
        freq[s, :len(c)] = [q[2] for q in c]
        drift[s, :len(c)] = [q[3] * weight for q in c]
    rows = np.arange(n_s)
    span = np.ptp(ts)
    nus = np.arange(1 / (40 * DAY), 1 / (0.15 * DAY), 0.05 / span)
    east, north = rdot @ E_EAST, rdot @ E_NORTH

    def one(nu):
        m = 2 * np.pi * nu * (ts - ts.mean())
        cols, dcols = [np.ones_like(m), east, north], [np.zeros_like(m)] * 3
        for k in range(1, harmonics + 1):
            cols += [np.cos(k * m), np.sin(k * m)]
            dcols += [-k * 2 * np.pi * nu * np.sin(k * m) * weight, k * 2 * np.pi * nu * np.cos(k * m) * weight]
        a = np.r_[np.array(cols).T, np.array(dcols).T]
        pick = np.zeros(n_s, int)
        for _ in range(rounds):
            y = np.r_[freq[rows, pick], drift[rows, pick]]
            res = y - a @ np.linalg.lstsq(a, y, rcond=None)[0]
            d2 = res[:n_s] ** 2 + res[n_s:] ** 2
            good = d2 <= np.sort(d2)[int(0.8 * n_s) - 1]
            both = np.r_[good, good]
            model = a @ np.linalg.lstsq(a[both], y[both], rcond=None)[0]
            pick = np.argmin((freq - model[:n_s, None]) ** 2 + (drift - model[n_s:, None]) ** 2, axis=1)
        y = np.r_[freq[rows, pick], drift[rows, pick]]
        return float(np.sum((y - a @ np.linalg.lstsq(a, y, rcond=None)[0]) ** 2)), pick

    scored = [one(nu) for nu in nus]
    score = np.array([q[0] for q in scored])
    picked = []
    for i in np.argsort(score):
        if all(abs(i - j) > 8 for j in picked):
            picked.append(int(i))
        if len(picked) == keep:
            break
    trial = {}
    for i in picked:
        for div in (1, 2, 3):
            nu = nus[i] * div
            if nu < nus[-1]:
                lo = max(np.searchsorted(nus, nu) - 3, 0)
                j = lo + int(np.argmin(score[lo:lo + 7]))
                trial[j] = (score[j], 1 / nus[j], scored[j][1])
    return sorted(trial.values(), key=lambda q: q[0])[:most]


def doppler_model(q, ts, rdot):
    """apparent frequency and drift at the session times for q = (f0, pb, x, e, om, t0, xi, eta)"""
    def freq(t):
        T = emission(t, q[1:6])
        return q[0] / (1 + roemer_rate(T, *q[1:6]))
    pos = 1 + (rdot @ E_EAST) * q[6] * ARCSEC + (rdot @ E_NORTH) * q[7] * ARCSEC
    return freq(ts) * pos, (freq(ts + 5.0) - freq(ts - 5.0)) / 10.0


def orbit_fit(q, ts, f, fd, rdot, weight, iters):
    q = q.copy()

    def res(qq):
        a, b = doppler_model(qq, ts, rdot)
        return np.r_[f - a, (fd - b) * weight]

    def jacobian_q(qq):
        jac = np.empty((2 * len(ts), 8))
        for j in range(8):
            d = np.zeros(8)
            d[j] = QSTEP[j]
            jac[:, j] = (res(qq - d) - res(qq + d)) / (2 * QSTEP[j])
        return jac

    r0 = res(q)
    c2, lam = r0 @ r0, 1e-3
    for _ in range(iters):
        jac = jacobian_q(q)
        scale = np.sqrt((jac ** 2).sum(0)) + 1e-300
        js = jac / scale
        trial = q + np.linalg.solve(js.T @ js + lam * np.eye(8), js.T @ r0) / scale
        trial[3] = min(max(trial[3], 1e-4), 0.9)
        r1 = res(trial)
        if np.isfinite(r1 @ r1) and r1 @ r1 < c2:
            q, r0, c2, lam = trial, r1, r1 @ r1, lam / 3
        else:
            lam *= 5
    jac = jacobian_q(q)
    return q, c2, np.linalg.pinv(jac.T @ jac) * c2 / max(len(r0) - 8, 1)


def orbits(ts, cands, rdot, weight=600.0, keep=4):
    """Elliptical orbits through one counting per session, best first: (misfit, q, covariance, countings)"""
    n_s = len(ts)
    take = lambda pick: (np.array([cands[s][pick[s]][2] for s in range(n_s)]), np.array([cands[s][pick[s]][3] for s in range(n_s)]))

    def distance(q, pick):
        mf, md = doppler_model(q, ts, rdot)
        f, fd = take(pick)
        return (f - mf) ** 2 + ((fd - md) * weight) ** 2, mf, md

    out = []
    for _, pb, pick in harmonic_scan(ts, cands, rdot, weight=weight, keep=keep):
        f, fd = take(pick)
        amp = (f.max() - f.min()) / 2 / f.mean()
        best = None
        for e0 in (0.05, 0.3):
            for om0 in np.arange(4) * np.pi / 2:
                for ph in np.arange(6) / 6:
                    q0 = np.array([f.mean(), pb, amp * pb / (2 * np.pi), e0, om0, ts.mean() + ph * pb, 0.0, 0.0])
                    q, c2, cov = orbit_fit(q0, ts, f, fd, rdot, weight, iters=8)
                    if best is None or c2 < best[1]:
                        best = (q, c2, cov)
        q = best[0]
        for _ in range(4):                               # fit without the two worst sessions, then recount all
            d2, _, _ = distance(q, pick)
            good = d2 <= np.sort(d2)[n_s - 3]
            q = orbit_fit(q, ts[good], f[good], fd[good], rdot[good], weight, iters=15)[0]
            _, mf, md = distance(q, pick)
            again = np.array([int(np.argmin([(c[2] - mf[s]) ** 2 + ((c[3] - md[s]) * weight) ** 2 for c in cands[s]]))
                              for s in range(n_s)])
            same = np.array_equal(again, pick)
            pick = again
            f, fd = take(pick)
            if same:
                break
        q, c2, cov = orbit_fit(q, ts, f, fd, rdot, weight, iters=30)
        out.append((c2, q, cov, pick))
    out.sort(key=lambda z: z[0])
    return out


def jacobian(t, r, p):
    T = emission(t + (r @ E_EAST) * p[3] * ARCSEC + (r @ E_NORTH) * p[4] * ARCSEC, p[5:])
    jac = np.empty((len(t), 10))
    jac[:, 0], jac[:, 1], jac[:, 2] = 1.0, T, T * T / 2
    jac[:, 3], jac[:, 4] = p[1] * (r @ E_EAST) * ARCSEC, p[1] * (r @ E_NORTH) * ARCSEC
    for j in range(5, 10):
        d = np.zeros(10)
        d[j] = STEPS[j]
        jac[:, j] = (phase(t, r, p + d) - phase(t, r, p - d)) / (2 * STEPS[j])
    return jac


def refit(p, t, r, sig, n, prior_mean, prior_icov, iters=3):
    """Gauss-Newton on the turn counts n with a Gaussian prior; returns p, chi-square of the data, covariance"""
    for it in range(iters + 1):
        jac = jacobian(t, r, p) / sig[:, None]
        res = (n - phase(t, r, p)) / sig
        normal = jac.T @ jac + prior_icov
        scale = np.sqrt(np.diag(normal))
        if it == iters:
            break
        rhs = jac.T @ res + prior_icov @ (prior_mean - p)
        p = p + np.linalg.lstsq(normal / np.outer(scale, scale), rhs / scale, rcond=None)[0] / scale
        p[7] = min(max(p[7], 1e-5), 0.95)
    cov = np.linalg.pinv(normal / np.outer(scale, scale)) / np.outer(scale, scale)
    return p, float(res @ res), cov


def lll(basis, delta=0.75):
    """Lenstra-Lenstra-Lovasz reduction of the columns of `basis`; returns the reduced basis and the unimodular U"""
    b = basis.astype(float).copy()
    n = b.shape[1]
    u = np.eye(n, dtype=object)

    def gram_schmidt():
        q, mu, norm = np.zeros_like(b), np.zeros((n, n)), np.zeros(n)
        for i in range(n):
            v = b[:, i].copy()
            for j in range(i):
                mu[i, j] = b[:, i] @ q[:, j] / norm[j]
                v -= mu[i, j] * q[:, j]
            q[:, i], norm[i] = v, v @ v
        return mu, norm

    mu, norm = gram_schmidt()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            c = np.round(mu[k, j])
            if c != 0:
                b[:, k] -= c * b[:, j]
                u[:, k] -= int(c) * u[:, j]
                mu, norm = gram_schmidt()
        if norm[k] >= (delta - mu[k, k - 1] ** 2) * norm[k - 1]:
            k += 1
        else:
            b[:, [k, k - 1]] = b[:, [k - 1, k]]
            u[:, [k, k - 1]] = u[:, [k - 1, k]]
            mu, norm = gram_schmidt()
            k = max(k - 1, 1)
    return b, u


def nearest(tri, c, margin, most, bound=None):
    """whole vectors z with |tri z - c|^2 within `margin` of the smallest value (tri upper triangular); with
    `bound` given, all z below that fixed bound (at most `most` of them)"""
    n = len(c)
    z = np.zeros(n)
    for k in range(n - 1, -1, -1):                       # nearest-plane point: a first bound
        z[k] = np.round((c[k] - tri[k, k + 1:] @ z[k + 1:]) / tri[k, k])
    fixed = bound is not None
    bound = [bound if fixed else float(np.sum((tri @ z - c) ** 2)) + margin + 1e-9]
    found, z = [], np.zeros(n)

    def descend(k, partial):
        if k < 0:
            found.append((partial, z.copy()))
            if not fixed and partial + margin < bound[0]:
                bound[0] = partial + margin + 1e-9
            return
        centre = (c[k] - tri[k, k + 1:] @ z[k + 1:]) / tri[k, k]
        half = np.sqrt(max(bound[0] - partial, 0.0)) / abs(tri[k, k])
        for v in sorted(range(int(np.ceil(centre - half)), int(np.floor(centre + half)) + 1), key=lambda w: abs(w - centre)):
            d = partial + (tri[k, k] * (v - centre)) ** 2
            if d <= bound[0] and len(found) < most:
                z[k] = v
                descend(k - 1, d)

    descend(n - 1, 0.0)
    found.sort(key=lambda q: q[0])
    return [q for q in found if q[0] <= found[0][0] + margin + 1e-9]


def connect(data, groups, counts, q, qcov, margin=150.0, top=40, most=50000, info=None, reduce=True, search=True,
            pin_position=False, no_f1=False):
    """turn counts between sessions as an integer least-squares problem.

    Around the solution of step 2 the phase is linear in the ten parameters; with the parameters eliminated the
    misfit is a quadratic form in the whole numbers N (one per session, the middle session fixed). The lattice is
    reduced and every N within `margin` of the best is listed; each is then refitted with the full model."""
    t, r, err = data["t"], data["r"], data["err"]
    p0 = np.array([0.0, q[0], 0.0, q[6], q[7], q[1], q[2], q[3], q[4], q[5]])
    sig = p0[1] * err
    wide = np.linalg.inv(np.diag([1e3, 1e-3, 1e-13, 120.0, 120.0, 1e3, 1.0, 0.5, 6.0, 1e5]) ** 2)
    idx = [1, 5, 6, 7, 8, 9, 3, 4]
    prior = wide.copy()
    prior[np.ix_(idx, idx)] += np.linalg.pinv(25 * qcov + np.diag((QSTEP * 1e-3) ** 2))
    if no_f1:                                                             # constant rotation frequency
        prior[2, 2] += 1e40
    if pin_position:                                                      # pulsar taken at the catalogue position
        p0[3:5] = 0.0
        prior[3, 3] += 1e12
        prior[4, 4] += 1e12
    mids = np.array([t[g].mean() for g in groups])
    anchor = int(np.argmin(np.abs(mids - np.median(mids))))
    others = [s for s in range(len(groups)) if s != anchor]
    base = phase(t, r, p0)
    jac = jacobian(t, r, p0)
    whole = np.zeros(len(t))
    design = np.zeros((len(t), len(others)))
    for s, g in enumerate(groups):
        whole[g] = counts[s] + np.round(base[g][0])
        if s != anchor:
            design[g, others.index(s)] = 1.0
    y = (base - whole) / sig
    scale = np.sqrt(np.sum((jac / sig[:, None]) ** 2, 0) + np.diag(prior))
    val, vec = np.linalg.eigh(prior / np.outer(scale, scale))
    root = (vec * np.sqrt(np.clip(val, 0.0, None))).T                 # root.T @ root = scaled prior
    a_p = np.r_[-(jac / sig[:, None]) / scale, root]
    a_n = np.r_[design / sig[:, None], np.zeros((10, len(others)))]
    y = np.r_[(base - whole) / sig, np.zeros(10)]
    q_p, r_p = np.linalg.qr(a_p)
    clean = lambda v: v - q_p @ (q_p.T @ v)                           # parameters eliminated
    q_n, upper = np.linalg.qr(clean(a_n))
    target = q_n.T @ clean(y)
    reduced, u = lll(upper) if reduce else (upper, np.eye(len(others), dtype=object))
    qmat, tri = np.linalg.qr(reduced)
    if search:
        best = nearest(tri, qmat.T @ target, 0.0, most)[0][0]
        listed = nearest(tri, qmat.T @ target, margin, most, bound=best + margin)
    else:                                                             # rounding one after another, no search
        c, z = qmat.T @ target, np.zeros(len(others))
        for k in range(len(others) - 1, -1, -1):
            z[k] = np.round((c[k] - tri[k, k + 1:] @ z[k + 1:]) / tri[k, k])
        listed = [(float(np.sum((tri @ z - c) ** 2)), z)]
    if info is not None:
        info.append(dict(listed=len(listed), order=others))
    out = []
    for value, z in listed[:top]:
        n_s = np.array([int(v) for v in (u @ np.array([int(v) for v in z], dtype=object))], float)
        start = p0 + np.linalg.solve(r_p, q_p.T @ (y - a_n @ n_s)) / scale
        counts_all = whole + design @ n_s
        p, c2, cov = refit(start, t, r, sig, counts_all, p0 if (pin_position or no_f1) else start,
                           prior if (pin_position or no_f1) else wide, iters=4)
        if np.isfinite(c2):
            out.append((c2, p, cov, tuple(n_s.astype(int))))
    out.sort(key=lambda z: z[0])
    return out


def solve(data, f_mid, keep=3, first_only=False, **kw):
    t = data["t"]
    groups = sessions(t)
    cands = [local(t[g], f_mid) for g in groups]
    if first_only:
        cands = [c[:1] for c in cands]
    ts = np.array([t[g].mean() for g in groups])
    rd = np.array([np.polyfit(t[g] - t[g].mean(), data["r"][g], 1)[0] for g in groups])
    found = orbits(ts, cands, rd)
    sols, lattice = [], []
    for c2, q, qcov, pick in found:
        if c2 > 25 * found[0][0] or len(lattice) >= keep or not (np.all(np.isfinite(q)) and np.all(np.isfinite(qcov)) and q[1] > 0):
            continue
        counts = [cands[s][pick[s]][1] for s in range(len(groups))]
        try:
            sols += connect(data, groups, counts, q, qcov, info=lattice, **kw)
        except np.linalg.LinAlgError:
            continue
    sols.sort(key=lambda z: z[0])
    info = dict(cands=cands, orbit=[(c[0], c[1][1] / DAY) for c in found], lattice=lattice, groups=groups,
                picks=[c[3] for c in found])
    return sols, info


def load(folder):
    rows = np.genfromtxt(Path(folder) / "toas.csv", delimiter=",", names=True)
    t = (rows["mjd_day"] - T_REF) * DAY + rows["sec_of_day"]
    order = np.argsort(t)
    r = np.c_[rows["x_ls"], rows["y_ls"], rows["z_ls"]]
    return dict(t=t[order], err=rows["err_us"][order] * 1e-6, r=r[order])


def results(p):
    p = np.array(p, float)
    if p[6] < 0:                                         # same orbit with x > 0
        p[6], p[8] = -p[6], p[8] + np.pi
    pb = p[5] / DAY
    t0 = p[9] / DAY
    return {"spin_frequency_hz": float(p[1]), "offset_east_arcsec": float(p[3]), "offset_north_arcsec": float(p[4]), "orbital_period_days": float(pb),
            "projected_semimajor_axis_ls": float(p[6]), "eccentricity": float(p[7]),
            "periastron_longitude_deg": float(np.degrees(p[8]) % 360.0),
            "periastron_epoch_mjd": float(T_REF + t0 - pb * np.floor(t0 / pb))}


def main():
    data = load(sys.argv[1])
    sols, info = solve(data, F_APPARENT)
    out = results(sols[0][1])
    out["spin_frequency_derivative_hz_per_s"] = float(sols[0][1][2])
    out["chi_square"], out["arrival_times"] = sols[0][0], int(len(data["t"]))
    out["chi_square_of_the_runner_up"] = min((s[0] for s in sols[1:] if s[3] != sols[0][3]), default=None)
    out["turn_counts_within_150_of_the_best"] = [q["listed"] for q in info["lattice"]]
    Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
    Path(sys.argv[2]).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out))


if __name__ == "__main__":
    main()
