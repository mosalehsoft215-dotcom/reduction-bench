"""Recover the Pearl length and signed vortex locations from a flux survey."""
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import j0, j1, roots_legendre


def load_approach(folder, head, *, stage=True):
    """Height, cold stage scale and the base-temperature Pearl length from the two approaches."""
    data = {"reference": [], "film": []}
    with (folder / "approaches.csv").open() as stream:
        for row in csv.DictReader(stream):
            data[row["surface"]].append([float(row["lift_um"]),
                float(row["screen_on_pH"]) - float(row["normal_background_pH"])])
    root, weight = roots_legendre(190)
    k, dk = (root + 1) * 2.7, weight * 2.7
    a, b = head["drive_radius_um"], head["pickup_radius_um"]
    spectral = -4e-7 * np.pi**2 * 1e6 * a*b * dk*j1(k*a)*j1(k*b)
    blocks = {}
    for name, rows in data.items():
        dz, observed = np.asarray(rows).T
        basis = np.column_stack((np.ones_like(dz), dz))
        q = np.linalg.qr(basis, mode="reduced")[0]
        offset = head["reference_screen_offset_um"] if name == "reference" else 0.
        blocks[name] = q, observed, dz, offset

    def residual(params):
        height, scale, lam = params
        rows = []
        for name, (q, observed, dz, offset) in blocks.items():
            transfer = spectral if name == "reference" else spectral/(1+k*lam)
            theory = np.exp(-2*(height + scale*dz[:,None] + offset)*k) @ transfer
            difference = observed - theory
            rows.append(difference - q @ (q.T @ difference))
        return np.concatenate(rows)

    lo, hi = np.asarray([head["height_bounds_um"], (.5, 1.5) if stage else (1., 1.0000001),
                         head["pearl_bounds_um"]]).T
    best = None
    for start in (.65, .85, 1.0, 1.2):
        x0 = np.clip(np.r_[(lo[0]+hi[0])/2, start, (lo[2]+hi[2])/2], lo, hi)
        fit = least_squares(residual, x0, bounds=(lo,hi), x_scale="jac", ftol=1e-12, xtol=1e-12, gtol=1e-12)
        if best is None or fit.cost < best.cost:
            best = fit
    return best.x


class ImageModel:
    def __init__(self, folder, head, height, lam, *, scale=1., unwrap=True, surroundings=True, finite_loop=True):
        with np.load(folder / "rasters.npz") as arrays:
            self.xy = arrays["xy_um"]
            flux = arrays["wrapped_phi0"]
            if unwrap:
                flux = np.unwrap(flux*2*np.pi, axis=-1)/(2*np.pi)
            self.observed = self.project(flux)
            z = scale*arrays["lift_um"] + height
            interval = arrays["sample_interval_s"]
        # The flux readout and the encoders are latched by separate trigger lines: a flux sample belongs to where
        # the loop was a fixed delay before the encoder latch, i.e. back along the row by velocity x delay.
        self.encoder = self.xy.copy()
        self.travel = np.diff(self.encoder, axis=2).mean(axis=2, keepdims=True) / interval[:, None, None, None]
        direction = np.sign(self.travel[..., 0].mean(axis=(1, 2)))
        self.opposed = bool(direction.min() != direction.max())
        self.delay = 0.
        self.noise = head["flux_sigma_phi0"]
        region = "search_bounds_um" if surroundings else "reported_bounds_um"
        self.bounds = np.asarray(head[region]).reshape(2,2)
        n, w = roots_legendre(220)
        self.k, self.dk = (n+1)*2.7, w*2.7
        b = head["pickup_radius_um"]
        self.radius = np.linspace(0., 120., 2401)
        self.shape = j0(self.radius[:,None]*self.k)
        self.pickup = b*j1(b*self.k) if finite_loop else b*b*self.k/2
        self.z = z
        self.set_pearl(lam)

    def set_delay(self, delay):
        self.delay = float(delay)
        self.xy = self.encoder - self.travel*delay

    def set_pearl(self, lam):
        spectral = self.dk*self.pickup/(1+self.k*lam)
        self.kernels = [CubicSpline(self.radius, self.shape @ (spectral*np.exp(-self.k*h))) for h in self.z]

    @staticmethod
    def project(a):
        t = np.linspace(-1., 1., a.shape[-1])
        return a - a.mean(axis=-1, keepdims=True) - (a @ t)[...,None]*t/(t@t)

    def columns(self, position):
        q = self.xy[None] - np.asarray(position)[:,None,None,None]
        r = np.linalg.norm(q, axis=-1)
        return self.project(np.stack([kernel(r[:, i]) for i,kernel in enumerate(self.kernels)], axis=1))

    def evaluate(self, params, signs, jac=False, free_delay=False):
        # A single common normalization absorbs first-order uncertainty
        # in the approach calibration during catalogue reconstruction.
        if free_delay:
            self.set_delay(params[-1])
            params = params[:-1]
        gain = params[-1]
        q = self.xy[None] - params[:-1].reshape(-1,2)[:,None,None,None]
        r = np.linalg.norm(q, axis=-1)
        fields = np.stack([kernel(r[:, i]) for i,kernel in enumerate(self.kernels)], axis=1)
        combined = self.project(np.einsum("n,nijk->ijk", signs, fields))
        if not jac:
            return (gain*combined - self.observed).ravel()
        slope = np.stack([kernel(r[:,i], 1) for i,kernel in enumerate(self.kernels)], axis=1)
        gradients = (-gain*signs[:,None,None,None,None]*q*slope[...,None]
                     / np.maximum(r[...,None], 1e-10))
        spatial = self.project(np.moveaxis(gradients,-1,1)).reshape(len(signs)*2,-1).T
        blocks = [spatial, combined.ravel()]
        if free_delay:      # moving every sample back along its travel is the same as moving the sources forward
            along = np.einsum("nsijk,sijk->sij", gradients, np.broadcast_to(self.travel, gradients.shape[1:]))
            blocks.append(self.project(along).ravel())
        return np.column_stack(blocks)

    def fit(self, positions, signs, spread=.04, effort=65, free_delay=False):
        free_delay = free_delay and self.opposed
        x0 = np.r_[np.ravel(positions),1.]
        low = np.r_[np.tile(self.bounds[:,0],len(signs)),1-spread]
        high = np.r_[np.tile(self.bounds[:,1],len(signs)),1+spread]
        if free_delay:
            x0, low, high = np.r_[x0, self.delay], np.r_[low, self.delay-.1], np.r_[high, self.delay+.1]
        fit = least_squares(self.evaluate,x0,jac=lambda p,s,free_delay:self.evaluate(p,s,True,free_delay),
                            args=(signs,),kwargs={"free_delay":free_delay},bounds=(low,high),
                            max_nfev=effort,ftol=1e-8,xtol=1e-8)
        if free_delay:
            self.set_delay(fit.x[-1])
            fit.x = fit.x[:-1]
        return fit

    def recover(self):
        x,y = [np.linspace(lo+.1,hi-.1, int((hi-lo)/2.4)+1) for lo,hi in self.bounds]
        grid = np.stack(np.meshgrid(x,y),axis=-1).reshape(-1,2)
        dictionary = np.concatenate([self.columns(chunk).reshape(len(chunk),-1).astype("float32")
                                     for chunk in np.array_split(grid, 16)])
        dictionary /= np.linalg.norm(dictionary,axis=1)[:,None]
        residual = -self.observed.ravel()
        locations = np.empty((0,2))
        signs = np.empty(0)
        for _ in range(11):
            scores = -(dictionary @ residual.astype("float32"))
            index = int(np.argmax(abs(scores)))
            if abs(scores[index]) < 10*self.noise:
                break
            signs = np.append(signs,np.sign(scores[index]))
            # Intermediate catalogues only seed the next pick; the last one is refined in full below.
            fit = self.fit(np.vstack((locations,grid[index])),signs,effort=18)
            locations, residual = fit.x[:-1].reshape(-1,2), fit.fun
        locations, signs = self.drop_cancelled_pairs(locations, signs, residual)
        if len(signs):
            final = self.fit(locations, signs)
            locations = final.x[:-1].reshape(-1,2)
            self.misfit = float(final.fun @ final.fun)
        else:
            self.misfit = float(self.observed.ravel() @ self.observed.ravel())
        return locations, signs

    def drop_cancelled_pairs(self, locations, signs, residual):
        # A fitted pair of almost coincident opposite signs can absorb
        # noise but leave no resolved vortex signal. Check the reduced fit.
        while len(signs) > 1:
            separation = np.linalg.norm(locations[:,None] - locations[None,:],axis=-1)
            possible = (signs[:,None] != signs[None,:]) & (separation < 1.5)
            if not possible.any():
                break
            pair = np.unravel_index(np.argmin(np.where(possible,separation,np.inf)),separation.shape)
            keep = np.ones(len(signs),dtype=bool)
            keep[list(pair)] = False
            if keep.any():
                fit = self.fit(locations[keep],signs[keep])
                trial, difference = fit.x[:-1].reshape(-1,2), fit.fun
            else:
                trial, difference = locations[keep], -self.observed.ravel()
            if difference@difference - residual@residual > 80*self.noise**2:
                break
            locations, signs, residual = trial, signs[keep], difference
        return locations, signs


def latch_delay(model, span=.25):
    """Trigger delay from rasters run in opposite directions: stationary vortices must line up once each raster's
    samples are moved back along their travel. Rasters at different lifts differ only by a symmetric blur, so the
    normalised cross-correlation of opposite-direction pairs peaks at the true delay."""
    direction = np.sign(model.travel[..., 0].mean(axis=(1, 2)))
    pairs = [(a, b) for a in range(len(direction)) for b in range(a+1, len(direction)) if direction[a] != direction[b]]
    if not pairs:
        return 0.                     # equal travel everywhere: the delay is a common shift the data cannot see
    rows, cols = model.observed.shape[1:]
    j, i = np.mgrid[:rows, :cols]
    basis = np.column_stack((i.ravel(), j.ravel(), np.ones(i.size)))
    affine = [np.linalg.lstsq(basis, e.reshape(-1, 2), rcond=None)[0] for e in model.encoder]
    velocity = model.travel.mean(axis=(1, 2, 3))
    lo, hi = model.encoder.min(axis=(0, 1, 2)) + 3., model.encoder.max(axis=(0, 1, 2)) - 3.
    gx, gy = np.meshgrid(np.arange(lo[0], hi[0], .5), np.arange(lo[1], hi[1], .5))
    points = np.column_stack((gx.ravel(), gy.ravel()))

    def sample(s, delay):
        a = affine[s]
        ij = np.linalg.solve(a[:2].T, (points + velocity[s]*delay - a[2]).T)
        return map_coordinates(model.observed[s], [ij[1], ij[0]], order=3, mode="constant", cval=np.nan)

    def score(delay):
        images = [sample(s, delay) for s in range(len(direction))]
        total = 0.
        for a, b in pairs:
            ok = np.isfinite(images[a]) & np.isfinite(images[b])
            u, w = images[a][ok] - images[a][ok].mean(), images[b][ok] - images[b][ok].mean()
            total += u @ w / np.sqrt((u @ u) * (w @ w))
        return -total

    trial = np.linspace(-span, span, 201)
    best = trial[int(np.argmin([score(d) for d in trial]))]
    step = trial[1] - trial[0]
    return float(minimize_scalar(score, bounds=(best-step, best+step), method="bounded",
                                 options={"xatol": 1e-5}).x)


def scan_pearl(model, coordinates, signs, guess, bounds, span=(.6, 1.8)):
    """Pearl length at the scan temperature: the flux quantum fixes each image's absolute size, so the
    three lifts settle the screening that the base-temperature approach cannot."""
    def misfit(log_lam):
        model.set_pearl(np.exp(log_lam))
        fit = model.fit(coordinates, signs, spread=1e-6)
        return fit.fun @ fit.fun
    lo, hi = np.log(bounds)
    grid = np.linspace(max(lo, np.log(guess*span[0])), min(hi, np.log(guess*span[1])), 5)
    values = [misfit(g) for g in grid]
    i = int(np.argmin(values))
    left, right = grid[max(i-1, 0)], grid[min(i+1, len(grid)-1)]
    best = minimize_scalar(misfit, bounds=(left, right), method="bounded", options={"xatol": 3e-4})
    model.set_pearl(np.exp(best.x))
    return float(np.exp(best.x))


def reduce(folder, *, stage=True, scan_temperature=True, trigger=True, budget=150., **image_options):
    start = time.monotonic()
    head = json.loads((folder/"head.json").read_text())
    height, scale, lam = load_approach(folder, head, stage=stage)
    model = ImageModel(folder,head,height,lam,scale=scale,**image_options)
    model.set_delay(latch_delay(model) if trigger else 0.)
    if trigger and model.opposed:
        # A first catalogue at the cross-correlation delay, then the delay refined with its coordinates.
        model.set_pearl(lam*1.12)
        coordinates, signs = model.recover()
        if len(signs):
            model.fit(coordinates, signs, free_delay=True)
        model.set_pearl(lam)
    if not scan_temperature:
        coordinates, signs = model.recover()
    else:
        # The approach ran at base temperature, colder than the rasters. Each trial screening length gets its
        # own catalogue; the scan length is the one whose catalogue leaves only noise (a count penalty
        # keeps spurious sources from paying for a wrong length).
        trials = []
        for factor in (1.12, 1.0, 1.26):            # most likely first: a slow host may stop early
            if trials and time.monotonic() - start > budget:     # a slow host keeps the best catalogue so far
                break
            model.set_pearl(lam*factor)
            found = model.recover()
            trials.append((model.misfit/model.noise**2 + 40*len(found[1]), lam*factor, found))
        score, lam, (coordinates, signs) = min(trials, key=lambda item: item[0])
        for _ in range(3):
            if time.monotonic() - start > 1.6*budget:
                model.set_pearl(lam)
                break
            lam = scan_pearl(model, coordinates, signs, lam, head["pearl_bounds_um"], span=(.9, 1.1))
            if len(signs):
                # The flux quantum fixes the image amplitude, so positions are refitted with unit gain.
                final = model.fit(coordinates, signs, spread=1e-6, free_delay=trigger)
                coordinates = final.x[:-1].reshape(-1, 2)
                score = final.fun @ final.fun / model.noise**2 + 40*len(signs)
            # A fresh search at the settled length replaces the catalogue only if it explains the data better.
            fresh, fresh_signs = model.recover()
            if model.misfit/model.noise**2 + 40*len(fresh_signs) >= score - 1:
                break
            coordinates, signs = fresh, fresh_signs
    xl,xh,yl,yh = head["reported_bounds_um"]
    vortices = [{"x_um":float(x),"y_um":float(y),"winding":int(sign)}
                for (x,y),sign in zip(coordinates,signs) if xl<=x<=xh and yl<=y<=yh]
    return {"pearl_um":float(lam),"zero_lift_um":float(height),"vortices":vortices}


if __name__ == "__main__":
    result=reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result,allow_nan=False)+"\n", encoding="utf8")
