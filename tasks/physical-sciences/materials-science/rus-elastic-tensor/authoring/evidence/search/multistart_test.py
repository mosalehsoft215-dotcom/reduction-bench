import sys, time
import numpy as np
from scipy.stats import qmc
sys.path.insert(0, "/b/solution")
import rus_invert as ri


def stable(c):
    """Born criteria for a 4/mmm tensor in the sample frame."""
    c11, c33, c44, c66, c12, c13 = c
    return (c44 > 0 and c66 > 0 and c11 > abs(c12)
            and c33 * (c11 + c12) > 2 * c13 ** 2)

TRUE_S = np.array([162.45, 139.6, 52.3, 53.75, 72.85, 82.7])
t0 = time.time()
dims, rho = ri.read_sample("/b/environment/data/sample.txt")
pa, _ = ri.extract_peaks("/b/environment/data/mount_A.dat")
pb, _ = ri.extract_peaks("/b/environment/data/mount_B.dat")
f_all = np.array([r["f"] for r in ri.merged_resonances(pa, pb)])
fm = f_all[:20]
m8 = ri.RusModel(dims, rho, 8)
print("setup %.1fs" % (time.time() - t0), flush=True)

def robust_score(c, model=m8, f=fm, cap=0.02):
    if not stable(c): return 99.0
    fp = model.spectrum(c, f.size + 12)
    tot, rows, cols = ri.assignment_cost(f, fp, cap)
    unmatched = np.sum(fp < f[-1]) - np.sum(fp[cols] < f[-1])
    return (tot + 0.3 * cap**2 * unmatched) / cap**2

def local_robust(c0, model=m8, f=fm, iters=16):
    c = np.array(c0, float)
    caps = (0.06, 0.04, 0.025, 0.015)
    for it in range(iters):
        cap = caps[min(it // 4, 3)]
        fp, J = model.spectrum(c, f.size + 12, jac=True)
        _, rows, cols = ri.assignment_cost(f, fp, cap)
        good = np.abs(np.log(f[rows] / fp[cols])) < cap
        if good.sum() < 8: break
        r = np.log(f[rows[good]]) - np.log(fp[cols[good]])
        a = J[cols[good]] * c            # d ln f / d ln C
        h = a.T @ a
        step = np.linalg.solve(h + 1e-3 * np.trace(h) / 6 * np.eye(6), a.T @ r)
        c = c * np.exp(np.clip(step, -0.25, 0.25))
        if not stable(c): break
    return c

lo = np.array([60, 60, 15, 15, 1, 1.0]); hi = np.array([320, 320, 130, 130, 160, 160.0])
for n_start in (128, 256):
    t = time.time()
    pts = qmc.scale(qmc.Sobol(6, scramble=True, seed=7).random(n_start), lo, hi)
    ends, scores = [], []
    for p in pts:
        if not stable(p):
            continue
        c = local_robust(p)
        ends.append(c); scores.append(robust_score(c))
    ends = np.array(ends); scores = np.array(scores)
    err = np.abs(ends / TRUE_S - 1).max(axis=1)
    right = err < 0.03
    order = np.argsort(scores)
    print(f"{n_start} starts ({len(ends)} stable) in {time.time()-t:.1f}s: "
          f"{right.sum()} reach the right basin; best score {scores[order[0]]:.3f} "
          f"(right basin? {right[order[0]]}); best wrong score {scores[~right].min():.3f}", flush=True)
    print("  right-basin scores:", np.round(np.sort(scores[right])[:8], 3))
    print("  lowest 8 scores & err:", [(round(scores[i],3), round(err[i],3)) for i in order[:8]])
