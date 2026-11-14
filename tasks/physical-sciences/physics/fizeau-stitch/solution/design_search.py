"""Adaptive follow-up policy search used by reconstruct_surface.py.

The score is a max-min eigenvalue over seven scenarios with a discrete pilot, a data-dependent branch split
and three ordered six-pose routes, so the landscape has many local optima. Each seeded run anneals the full
policy jointly (pilot, routes, units, order; the branch split is re-chosen at every step among the balanced
splits that keep the 0.008 threshold margin), then applies an exhaustive best-improvement polish, a short
re-anneal and a second polish. Seeds are run two at a time until one clears the required score.
"""
import json, os
import numpy as np

_A = {}
PATS = [(a, b, 7 - a - b) for a in (1, 2, 3) for b in (1, 2, 3) if 1 <= 7 - a - b <= 3]
MARGIN = 0.0082   # verifier needs >= 0.008
EYE = 1e-6 * np.eye(8)


def init(infos, drift, diag, angles, dxs, dys, alpha, prior):
    _A.update(M=infos, G=drift, D=diag, ANG=angles, DX=dxs, DY=dys, AL=alpha, PR=prior, NC=diag.shape[1])


def seg(js, us):
    A = _A; js = np.asarray(js); us = np.asarray(us, float)
    pa = np.r_[0., A['ANG'][js[:-1]]]; px = np.r_[0., A['DX'][js[:-1]]]; py = np.r_[0., A['DY'][js[:-1]]]
    mv = np.abs(((A['ANG'][js] - pa + 180) % 360) - 180) / 90 + 0.75 * np.hypot(A['DX'][js] - px, A['DY'][js] - py) / 0.070710678
    eu = us[None, :] * (0.38 + 0.62 * np.exp(-A['AL'][:, None] * mv[None, :]))
    tm = np.cumsum(us) - 0.5 * us
    F = np.einsum('sk,skab->sab', eu, A['M'][:, js]); g = np.einsum('sk,ska->sa', eu * tm, A['G'][:, js])
    h = A['PR'] + (eu * tm ** 2).sum(1)
    return F - g[:, :, None] * g[:, None, :] / h[:, None, None]


class State:
    def __init__(s, pil, routes):
        s.pil = list(pil); s.routes = [[list(j), list(u)] for j, u in routes]; s.full()

    def full(s):
        s.P = seg(s.pil, [2, 2, 2]); s.sig = _A['D'][:, s.pil].mean(1)
        s.E = np.array([np.linalg.eigvalsh(EYE + s.P + seg(*r))[:, 0] for r in s.routes])

    def route_eigs(s, js, us):
        return np.linalg.eigvalsh(EYE + s.P + seg(js, us))[:, 0]

    def score(s, E=None, sig=None):
        E = s.E if E is None else E; sig = s.sig if sig is None else sig
        o = np.argsort(sig); ss = sig[o]; best = (-1.0, None)
        for a, b, c in PATS:
            if min(ss[a] - ss[a - 1], ss[a + b] - ss[a + b - 1]) / 2 < MARGIN: continue
            v = min(E[0, o[:a]].min(), E[1, o[a:a + b]].min(), E[2, o[a + b:]].min())
            if v > best[0]: best = (v, (a, b, c))
        return best


def _random_route(rng, excl):
    js = [int(j) for j in rng.choice([j for j in range(_A['NC']) if j not in excl], 6, replace=False)]; us = [3] * 6
    for _ in range(20):
        a, c = rng.integers(6, size=2)
        if a != c and us[a] > 1 and us[c] < 6: us[a] -= 1; us[c] += 1
    return js, us


def anneal(seed, iters, T0=3e-3, T1=2e-5, st=None):
    rng = np.random.default_rng(seed); NC = _A['NC']
    if st is None:
        pil = [int(j) for j in rng.choice(NC, 3, replace=False)]
        st = State(pil, [_random_route(rng, set(pil)) for _ in range(3)])
    cur = st.score()[0]; best = (cur, json.dumps([st.pil, st.routes]))
    for it in range(iters):
        T = T0 * (T1 / T0) ** (it / iters); r = rng.random()
        if r < 0.12:
            pil = st.pil.copy()
            if rng.random() < 0.3:
                i, k = rng.choice(3, 2, replace=False); pil[i], pil[k] = pil[k], pil[i]
            else:
                used = set(pil) | {j for rr in st.routes for j in rr[0]}
                pil[rng.integers(3)] = int(rng.choice([j for j in range(NC) if j not in used]))
            P = seg(pil, [2, 2, 2]); sig = _A['D'][:, pil].mean(1)
            E = np.array([np.linalg.eigvalsh(EYE + P + seg(*rr))[:, 0] for rr in st.routes]); v = st.score(E, sig)[0]
            if v >= cur or rng.random() < np.exp((v - cur) / T): st.pil, st.P, st.sig, st.E, cur = pil, P, sig, E, v
        else:
            b = int(rng.integers(3)); js, us = st.routes[b][0].copy(), st.routes[b][1].copy()
            if r < 0.55:
                excl = set(js) | set(st.pil); js[rng.integers(6)] = int(rng.choice([j for j in range(NC) if j not in excl]))
            elif r < 0.72:
                i, k = rng.choice(6, 2, replace=False); js[i], js[k] = js[k], js[i]
            else:
                a, c = rng.choice(6, 2, replace=False)
                if us[a] <= 1 or us[c] >= 6: continue
                us[a] -= 1; us[c] += 1
            E = st.E.copy(); E[b] = st.route_eigs(js, us); v = st.score(E)[0]
            if v >= cur or rng.random() < np.exp((v - cur) / T): st.routes[b] = [js, us]; st.E = E; cur = v
        if cur > best[0]: best = (cur, json.dumps([st.pil, st.routes]))
    pil, routes = json.loads(best[1]); return State(pil, routes)


def polish(st):
    NC = _A['NC']; cur = st.score()[0]
    while True:
        bm = (cur, None)
        for b in range(3):
            js, us = st.routes[b]; excl = set(js) | set(st.pil)
            trials = [(js[:i] + [j] + js[i + 1:], us) for i in range(6) for j in range(NC) if j not in excl]
            for i in range(6):
                for k in range(i + 1, 6):
                    nj = js.copy(); nj[i], nj[k] = nj[k], nj[i]; trials.append((nj, us))
            for a in range(6):
                for c in range(6):
                    if a != c and us[a] > 1 and us[c] < 6:
                        nu = us.copy(); nu[a] -= 1; nu[c] += 1; trials.append((js, nu))
            for nj, nu in trials:
                E = st.E.copy(); E[b] = st.route_eigs(nj, nu); v = st.score(E)[0]
                if v > bm[0] + 1e-12: bm = (v, ('r', b, nj, nu))
        used = {j for rr in st.routes for j in rr[0]}
        for i in range(3):
            for j in range(NC):
                if j in used or j in st.pil: continue
                pil = st.pil.copy(); pil[i] = j; v = State(pil, st.routes).score()[0]
                if v > bm[0] + 1e-12: bm = (v, ('p', pil))
        if bm[1] is None: return st, cur
        m = bm[1]
        if m[0] == 'r': st.routes[m[1]] = [m[2], m[3]]; st.full()
        else: st = State(m[1], st.routes)
        cur = st.score()[0]


def run_seed(seed, iters=300000):
    st = anneal(seed, iters); st, _ = polish(st)
    st = anneal(seed + 1000, iters // 3, T0=6e-4, st=st); st, v = polish(st)
    return v, st.pil, st.routes


def search_policy(cand_ids, gate, max_seeds=16, workers=2):
    best = None
    seeds = list(range(1, max_seeds + 1))
    for k in range(0, len(seeds), workers):
        batch = seeds[k:k + workers]
        if workers > 1 and hasattr(os, 'fork'):
            import multiprocessing as mp
            with mp.get_context('fork').Pool(len(batch)) as pool: res = pool.map(run_seed, batch)
        else:
            res = [run_seed(s) for s in batch]
        for r in res:
            if best is None or r[0] > best[0]: best = r
        print(f'design search: seeds {batch} -> best {best[0]:.6f}', flush=True)
        if best[0] >= gate + 5e-5: break
    v, pil, routes = best; st = State(pil, routes); val, (a, b, c) = st.score(); ss = np.sort(st.sig)
    names = ['low', 'mid', 'high']
    policy = {'pilot_plan': [{'candidate_id': int(cand_ids[j]), 'exposure_units': 2} for j in pil],
              'branch_thresholds': [float((ss[a - 1] + ss[a]) / 2), float((ss[a + b - 1] + ss[a + b]) / 2)],
              'followup_plans': {names[i]: [{'candidate_id': int(cand_ids[j]), 'exposure_units': int(u)}
                                            for j, u in zip(*routes[i])] for i in range(3)}}
    return policy, float(val)
