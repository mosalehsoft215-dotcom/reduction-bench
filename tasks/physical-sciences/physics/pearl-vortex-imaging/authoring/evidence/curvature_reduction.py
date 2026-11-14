"""Independent estimator based on curvature rasters and composite Simpson quadrature."""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import j0, j1

K = np.linspace(0, 5.8, 581)
W = np.ones(K.size)
W[1:-1:2], W[2:-1:2] = 4, 2
W *= (K[1] - K[0]) / 3


def calibrated_geometry(folder, settings):
    with (folder / "approaches.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    a,b = settings["drive_radius_um"],settings["pickup_radius_um"]
    spectral = -4*np.pi**2*.1*a*b*j1(K*a)*j1(K*b)*W

    def sweep(surface):
        raw = [r for r in rows if r["surface"] == surface]
        lifts = np.asarray([float(r["lift_um"]) for r in raw])
        values = np.asarray([float(r["screen_on_pH"]) - float(r["normal_background_pH"]) for r in raw])
        design = np.array([np.ones_like(lifts), lifts]).T
        return values, design, lifts

    def score(block, depth, pearl, scale=1.):
        values, design, lifts = block
        gap = values - np.exp(-2*scale*lifts[:,None]*K) @ (spectral * np.exp(-2*depth*K)/(1+K*pearl))
        gap -= design @ np.linalg.lstsq(design, gap, rcond=None)[0]
        return gap@ gap

    known, film = sweep("reference"), sweep("film")
    dz = settings["reference_screen_offset_um"]
    # The perfect screen alone fixes the height and the cold stage's displacement per commanded micrometre.
    def height_at(scale):
        return minimize_scalar(lambda z:score(known,z+dz,0,scale),bounds=settings["height_bounds_um"],
                               method="bounded",options={"xatol":1e-9})
    grid=np.linspace(.5,1.5,41)
    coarse=grid[int(np.argmin([height_at(s).fun for s in grid]))]
    scale=minimize_scalar(lambda s:height_at(s).fun,bounds=(coarse-.03,coarse+.03),
                          method="bounded",options={"xatol":1e-7}).x
    h=height_at(scale).x
    lam = minimize_scalar(lambda x:score(film,h,x,scale),bounds=settings["pearl_bounds_um"],
                          method="bounded",options={"xatol":1e-8}).x
    return h,scale,lam


def travel(folder):
    with np.load(folder/"rasters.npz",allow_pickle=False) as a:
        return np.diff(a["xy_um"],axis=2).mean(axis=2,keepdims=True)/a["sample_interval_s"][:,None,None,None]


def coordinates(folder, settings, height, pearl, scale, seed=None, delay=0., only=None):
    with np.load(folder/"rasters.npz",allow_pickle=False) as a:
        xy, lifts, wrapped = a["xy_um"], a["lift_um"], a["wrapped_phi0"]
    xy = xy - travel(folder)*delay      # flux latched a fixed delay before the encoders: back along the travel
    if only is not None:
        xy, lifts, wrapped = xy[only:only+1], lifts[only:only+1], wrapped[only:only+1]
    steps = np.diff(wrapped,axis=-1)
    increments = steps - np.rint(steps)
    target = np.diff(increments,axis=-1).ravel()
    r = np.arange(0.,121.,.07)
    b = settings["pickup_radius_um"]
    basis = j0(r[:,None]*K)
    transform = b*j1(b*K)*W/(1+pearl*K)
    tables = basis @ (transform[:,None]*np.exp(-K[:,None]*(height+scale*lifts)))
    lo_x,hi_x,lo_y,hi_y=settings["search_bounds_um"]
    px,py=np.meshgrid(np.arange(lo_x+.2,hi_x-.1,2.65),np.arange(lo_y+.2,hi_y-.1,2.65))
    candidates=np.column_stack((px.ravel(),py.ravel()))

    def dictionary(positions):
        rads = np.sqrt(np.sum((positions[:,None,None,None,:] - xy[None])**2,axis=-1))
        images = np.stack([np.interp(rads[:,i],r,tables[:,i]) for i in range(len(lifts))],axis=1)
        return np.diff(images,n=2,axis=-1).reshape(len(positions),-1)

    if seed is not None:      # refit a fixed catalogue under this screening and return its misfit
        charges=seed[1]
        fit=least_squares(lambda p:charges@dictionary(p.reshape(-1,2))-target,seed[0].ravel(),
                          diff_step=1e-5,max_nfev=20,ftol=1e-9)
        return fit.fun@fit.fun
    grid_columns=np.concatenate([dictionary(chunk) for chunk in np.array_split(candidates,20)])
    norms=np.linalg.norm(grid_columns,axis=1)
    model=np.empty((0,2))
    charges=np.empty(0)
    remaining=target.copy()
    sigma=settings["flux_sigma_phi0"]
    for _ in range(11):
        scores=grid_columns@remaining/norms
        ind=int(np.argmax(abs(scores)))
        if abs(scores[ind]) < 10*np.sqrt(6)*sigma:
            break
        charges=np.append(charges,np.sign(scores[ind]))
        initial=np.vstack((model,candidates[ind])).ravel()
        low=np.tile([lo_x,lo_y],len(charges))
        high=np.tile([hi_x,hi_y],len(charges))
        fit=least_squares(lambda p:charges@dictionary(p.reshape(-1,2))-target,initial,
                          bounds=(low,high),diff_step=1e-5,max_nfev=70,ftol=1e-8)
        model,remaining=fit.x.reshape(-1,2),-fit.fun
    return model,charges


def trigger_delay(folder, settings, h, lam, scale):
    """Separate catalogues of the first two rasters: when they ran in opposite directions, each matched vortex
    appears displaced forward along its raster's travel by velocity x delay."""
    v=travel(folder).mean(axis=(1,2))
    if v[0]@v[1] > 0:
        return 0.
    (p0,s0),(p1,s1)=[coordinates(folder,settings,h,lam,scale,only=i) for i in (0,1)]
    gap=v[0]-v[1]
    votes=[]
    for p,s in zip(p0,s0):
        d=np.linalg.norm(p1-p,axis=1)+1e9*(s1!=s)
        if len(d) and d.min()<6:
            votes.append((p-p1[int(np.argmin(d))])@gap/(gap@gap))
    delay=float(np.median(votes)) if votes else 0.
    positions,signs=coordinates(folder,settings,h,lam,scale,delay=delay)
    misfit=lambda t:coordinates(folder,settings,h,lam,scale,(positions,signs),delay=t)
    return float(minimize_scalar(misfit,bounds=(delay-.02,delay+.02),method="bounded",options={"xatol":2e-4}).x)


def reduce(folder):
    settings=json.loads((folder/"head.json").read_text())
    h,scale,lam=calibrated_geometry(folder,settings)
    delay=trigger_delay(folder,settings,h,lam,scale)
    positions,signs=coordinates(folder,settings,h,lam,scale,delay=delay)
    for _ in range(4):
        # The approach ran at base temperature; the three lifts of the rasters fix the scan-temperature length.
        lo,hi=np.log(settings["pearl_bounds_um"])
        misfit=lambda u:coordinates(folder,settings,h,np.exp(u),scale,(positions,signs),delay=delay)
        trial=np.linspace(max(lo,np.log(lam*.7)),min(hi,np.log(lam*1.6)),7)
        i=int(np.argmin([misfit(u) for u in trial]))
        u=minimize_scalar(misfit,bounds=(trial[max(i-1,0)],trial[min(i+1,6)]),method="bounded",
                          options={"xatol":3e-4}).x
        lam=float(np.exp(u))
        fresh,fresh_signs=coordinates(folder,settings,h,lam,scale,delay=delay)
        stable=len(fresh_signs)==len(signs)
        positions,signs=fresh,fresh_signs
        if stable:
            break
    xmin,xmax,ymin,ymax=settings["reported_bounds_um"]
    answer=[{"x_um":float(x),"y_um":float(y),"winding":int(s)} for (x,y),s in zip(positions,signs)
            if xmin<=x<=xmax and ymin<=y<=ymax]
    return {"pearl_um":float(lam),"zero_lift_um":float(h),"vortices":answer}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(Path(sys.argv[1])),allow_nan=False),encoding="utf8")
