"""Timing model of a binary pulsar (generator and reference).

Times are seconds of barycentric dynamical time counted from the reference epoch (MJD T_REF). A row of the data
gives the arrival time t_nom at the Solar-System barycentre computed for the catalogue position of the cluster,
and the vector r from the barycentre to the telescope in light-seconds (equatorial axes). With the pulsar xi
arcsec east and eta arcsec north of that position the arrival time at the barycentre is
    t = t_nom + r . (xi e_east + eta e_north).
The pulse was emitted at pulsar time T with t = T + D(T), D the Roemer delay of an elliptical orbit,
    D(T) = x [(cos E - e) sin w + sqrt(1 - e^2) sin E cos w],   E - e sin E = 2 pi (T - T0) / Pb,
and pulses leave when the rotational phase phi0 + F0 T + F1 T^2 / 2 is a whole number.
"""
import numpy as np

DAY = 86400.0
T_REF = 60500
ARCSEC = np.pi / 648000
RA0, DEC0 = np.radians(6.02), np.radians(-72.08)
E_EAST = np.array([-np.sin(RA0), np.cos(RA0), 0.0])
E_NORTH = np.array([-np.sin(DEC0) * np.cos(RA0), -np.sin(DEC0) * np.sin(RA0), np.cos(DEC0)])
NAMES = ("phi0", "f0", "f1", "xi", "eta", "pb", "x", "ecc", "om", "t0")


def anomaly(m, e):
    big = m + e * np.sin(m)
    for _ in range(10):
        big = big - (big - e * np.sin(big) - m) / (1 - e * np.cos(big))
    return big


def roemer(T, pb, x, e, om, t0):
    big = anomaly(2 * np.pi * (T - t0) / pb, e)
    return x * ((np.cos(big) - e) * np.sin(om) + np.sqrt(1 - e * e) * np.sin(big) * np.cos(om))


def roemer_rate(T, pb, x, e, om, t0):
    big = anomaly(2 * np.pi * (T - t0) / pb, e)
    return x * 2 * np.pi / pb * (-np.sin(big) * np.sin(om) + np.sqrt(1 - e * e) * np.cos(big) * np.cos(om)) / (1 - e * np.cos(big))


def emission(t, orbit):
    T = t - roemer(t, *orbit)
    for _ in range(5):
        T = t - roemer(T, *orbit)
    return T


def barycentric(t_nom, r, xi, eta):
    return t_nom + (r @ E_EAST) * xi * ARCSEC + (r @ E_NORTH) * eta * ARCSEC


def phase(t_nom, r, p):
    phi0, f0, f1, xi, eta = p[:5]
    T = emission(barycentric(t_nom, r, xi, eta), p[5:])
    return phi0 + f0 * T + 0.5 * f1 * T * T
