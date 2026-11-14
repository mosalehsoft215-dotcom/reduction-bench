"""The seven creep programmes (two shipped, five graded) and the writer that turns a simulation into the files the
analyst receives.  Run:  python3 programmes.py <output root> [seed offset]"""
import csv
import datetime as dt
import json
import os
import sys

import numpy as np

from creepmodel import R_GAS, Furnace, Slices, run_test, steady_rate

G_STD = 9.80665

DRAWINGS = {
    "CR6-30": dict(d=6.0, Lp=30.0, rho=6.0, D=10.0, span=41.0),
    "CR10-60": dict(d=10.0, Lp=60.0, rho=3.0, D=14.0, span=67.0),
    "CR8-40": dict(d=8.0, Lp=40.0, rho=4.0, D=12.0, span=50.0),
    "CR6-25": dict(d=6.0, Lp=25.0, rho=8.0, D=10.0, span=38.5),
}

# expansion: mean coefficient 20 C -> T (1/K); modulus in GPa
MATERIALS = {
    "800HT": dict(alpha=([20, 300, 500, 700, 800, 900], [14.4e-6, 15.9e-6, 16.6e-6, 17.3e-6, 17.7e-6, 18.0e-6]),
                  E=([20, 500, 700, 800, 900], [196, 170, 158, 151, 143]), n=4.73, Q=316.8e3, beta=2.6, ep=0.0024,
                  ref=(1.0e-8, 60.0, 725.0)),
    "2.25Cr-1Mo": dict(alpha=([20, 300, 400, 500, 600], [12.4e-6, 13.4e-6, 13.8e-6, 14.2e-6, 14.5e-6]),
                       E=([20, 400, 500, 600], [210, 185, 177, 166]), n=6.42, Q=356.3e3, beta=3.1, ep=0.0019,
                       ref=(1.0e-8, 160.0, 565.0)),
    "316H": dict(alpha=([20, 300, 500, 600, 700], [16.0e-6, 17.3e-6, 17.9e-6, 18.2e-6, 18.5e-6]),
                 E=([20, 500, 600, 700], [195, 165, 157, 149]), n=7.18, Q=411.6e3, beta=2.2, ep=0.0026,
                 ref=(1.0e-8, 175.0, 625.0)),
    "617": dict(alpha=([20, 500, 700, 800, 900, 1000], [11.6e-6, 13.6e-6, 14.4e-6, 14.9e-6, 15.4e-6, 15.9e-6]),
                E=([20, 700, 800, 900, 1000], [211, 171, 165, 157, 149]), n=4.91, Q=327.4e3, beta=2.9, ep=0.0022,
                ref=(1.0e-8, 38.0, 900.0)),
    "X10CrMoVNb9-1": dict(alpha=([20, 400, 500, 600, 700], [10.9e-6, 12.3e-6, 12.6e-6, 12.9e-6, 13.1e-6]),
                          E=([20, 500, 600, 700], [218, 187, 173, 158]), n=9.47, Q=503.2e3, beta=2.4, ep=0.0017,
                          ref=(1.0e-8, 165.0, 600.0)),
    "2618-T61": dict(alpha=([20, 100, 150, 200, 250], [21.8e-6, 22.3e-6, 22.7e-6, 23.1e-6, 23.5e-6]),
                     E=([20, 150, 200, 250], [74, 69, 66, 62]), n=5.23, Q=158.9e3, beta=3.4, ep=0.0025,
                     ref=(1.0e-8, 195.0, 200.0)),
    "Cu-OF": dict(alpha=([20, 200, 250, 300, 350], [16.7e-6, 17.3e-6, 17.5e-6, 17.7e-6, 17.9e-6]),
                  E=([20, 200, 300, 350], [125, 114, 107, 103]), n=4.36, Q=193.7e3, beta=3.0, ep=0.0027,
                  ref=(1.0e-8, 52.0, 300.0)),
}


def material(name):
    m = dict(MATERIALS[name])
    rate, s, T = m["ref"]
    m["lnA"] = np.log(rate) - m["n"] * np.log(s) + m["Q"] / (R_GAS * (T + 273.15))
    return m


def stress_for(m, rate, T_C):
    return float(np.exp((np.log(rate) - m["lnA"] + m["Q"] / (R_GAS * (T_C + 273.15))) / m["n"]))


# Each programme: material, design point, machines, furnaces, thermocouple certificates, the test list and the log.
# A test: id, machine, furnace, drawing, set temperature, target steady rate (1/s) used to pick the pan mass,
# and optional events.
def P(name, mat, design, machines, furnaces, tcs, tests, start, seed=0, notes=(), spare_tcs=(), true_ratio=None):
    """notes: (test, hours after that test ended, text) for log entries not tied to a simulated event."""
    return dict(name=name, mat=mat, design=design, machines=machines, furnaces=furnaces, tcs=tcs, tests=tests,
                start=start, seed=seed, notes=list(notes), spare_tcs=tuple(spare_tcs),
                true_ratio=dict(true_ratio or {}))


def tc_cert(base, slope):
    temps = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    return dict(reading_C=temps, correction_K=[round(base + slope * (t - 600) / 100, 2) for t in temps])


def T(tid, machine, furnace, drawing, set_C, rate, tcs, **ev):
    return dict(id=tid, machine=machine, furnace=furnace, drawing=drawing, set_C=set_C, rate=rate, tcs=tcs, ev=ev)


FURN = {
    "F1": Furnace(1.2, 0.010, 0.5, 0.004),
    "F2": Furnace(1.9, 0.016, -0.7, -0.006),
    "F3": Furnace(2.8, 0.022, 0.9, 0.008),
    "F4": Furnace(1.5, 0.020, 0.4, 0.003),
    "F5": Furnace(3.3, 0.025, -1.1, -0.004),
    "F6": Furnace(2.2, 0.030, 0.6, 0.005),
}

TCS = {
    "K-104": tc_cert(0.6, 0.12), "K-107": tc_cert(-0.4, 0.08), "K-111": tc_cert(0.9, -0.10),
    "K-118": tc_cert(-0.7, 0.15), "K-122": tc_cert(0.3, 0.05), "K-125": tc_cert(-1.1, 0.09),
    "K-131": tc_cert(1.2, -0.06), "K-136": tc_cert(-0.2, 0.11), "K-140": tc_cert(0.5, -0.12),
    "K-149": tc_cert(1.4, -0.13),
}

TRIO_A, TRIO_B, TRIO_C = ("K-104", "K-107", "K-111"), ("K-118", "K-122", "K-125"), ("K-131", "K-136", "K-140")


def programmes():
    out = {}
    out["P1"] = P("P1", "800HT", (712.0, 34.0), {"M2": 10.0, "M4": 10.0}, ["F1", "F2"], TRIO_A + TRIO_B, [
        T("T01", "M2", "F1", "CR6-30", 700, 2.6e-8, TRIO_A),
        T("T02", "M4", "F2", "CR6-30", 700, 9.0e-9, TRIO_B),
        T("T03", "M2", "F1", "CR6-30", 700, 3.4e-9, TRIO_A),
        T("T04", "M4", "F2", "CR6-30", 725, 2.8e-8, TRIO_B, excursion=(0.12, 0.17, 7.0)),
        T("T05", "M2", "F1", "CR6-30", 725, 9.5e-9, TRIO_A),
        T("T06", "M4", "F2", "CR6-30", 725, 3.6e-9, TRIO_B),
        T("T07", "M2", "F1", "CR6-30", 750, 2.5e-8, TRIO_A),
        T("T08", "M4", "F2", "CR6-30", 750, 8.8e-9, TRIO_B),
        T("T09", "M2", "F1", "CR6-30", 750, 3.3e-9, TRIO_A),
    ], dt.datetime(2024, 2, 5, 9, 30), seed=11)
    out["P2"] = P("P2", "2.25Cr-1Mo", (548.0, 118.0), {"M1": 10.0, "M6": 20.0}, ["F3", "F4"], TRIO_A + TRIO_C, [
        T("T01", "M6", "F3", "CR10-60", 540, 3.3e-9, TRIO_C),
        T("T02", "M1", "F4", "CR6-30", 540, 9.4e-9, TRIO_A),
        T("T03", "M1", "F4", "CR6-30", 540, 2.7e-8, TRIO_A),
        T("T04", "M6", "F3", "CR10-60", 565, 3.5e-9, TRIO_C, slip=(0.55, -0.0123)),
        T("T05", "M1", "F4", "CR6-30", 565, 9.1e-9, TRIO_A),
        T("T06", "M1", "F4", "CR6-30", 565, 2.9e-8, TRIO_A),
        T("T07", "M6", "F3", "CR10-60", 590, 3.2e-9, TRIO_C),
        T("T08", "M1", "F4", "CR6-30", 590, 9.6e-9, TRIO_A),
        T("T09", "M1", "F4", "CR6-30", 590, 2.6e-8, TRIO_A),
    ], dt.datetime(2023, 9, 18, 10, 0), seed=12)
    # h1: middle thermocouple of the second trio drifted (log only), hot curved furnace
    out["h1"] = P("h1", "316H", (608.0, 132.0), {"M3": 10.0, "M5": 10.0}, ["F5", "F1"], TRIO_A + TRIO_B, [
        T("T01", "M3", "F5", "CR6-25", 600, 2.7e-8, TRIO_B),
        T("T02", "M5", "F1", "CR6-25", 600, 8.9e-9, TRIO_A),
        T("T03", "M3", "F5", "CR6-25", 600, 3.1e-9, TRIO_B),
        T("T04", "M5", "F1", "CR6-25", 625, 2.6e-8, TRIO_A),
        T("T05", "M3", "F5", "CR6-25", 625, 9.3e-9, TRIO_B, drift={"K-122": 2.4}),
        T("T06", "M5", "F1", "CR6-25", 625, 3.4e-9, TRIO_A),
        T("T07", "M3", "F5", "CR6-25", 650, 2.4e-8, TRIO_B, drift={"K-122": 2.4}),
        T("T08", "M5", "F1", "CR6-25", 650, 9.8e-9, TRIO_A),
        T("T09", "M3", "F5", "CR6-25", 650, 3.2e-9, TRIO_B, drift={"K-122": 2.4}),
    ], dt.datetime(2024, 6, 3, 8, 45), seed=21, notes=[
        ("T09", 30, "K-122 checked against the reference thermometer at 625 C: it reads 2.4 K high. At its previous "
                    "check, before T05, it was in order; take the error as present in T05, T07 and T09."),
    ])
    # h2: mains failures inside the late part of two tests, a knocked extensometer, both drawings
    out["h2"] = P("h2", "617", (868.0, 26.4), {"M1": 10.0, "M6": 20.0}, ["F2", "F6"], TRIO_A + TRIO_C, [
        T("T01", "M6", "F6", "CR10-60", 850, 3.6e-9, TRIO_C),
        T("T02", "M1", "F2", "CR6-30", 850, 9.2e-9, TRIO_A, outage=(0.71, 9.5)),
        T("T03", "M1", "F2", "CR6-30", 850, 2.7e-8, TRIO_A),
        T("T04", "M6", "F6", "CR10-60", 900, 3.3e-9, TRIO_C, slip=(0.81, 0.0187)),
        T("T05", "M1", "F2", "CR6-30", 900, 9.6e-9, TRIO_A),
        T("T06", "M1", "F2", "CR6-30", 900, 2.5e-8, TRIO_A, outage=(0.78, 6.2)),
        T("T07", "M6", "F6", "CR10-60", 950, 3.4e-9, TRIO_C),
        T("T08", "M1", "F2", "CR6-30", 950, 9.0e-9, TRIO_A),
        T("T09", "M1", "F2", "CR6-30", 950, 2.8e-8, TRIO_A),
    ], dt.datetime(2024, 1, 15, 9, 0), seed=22)
    # h3: martensitic steel, long tests on the 10 mm drawing, one test stopped in primary creep
    out["h3"] = P("h3", "X10CrMoVNb9-1", (583.0, 147.0), {"M3": 10.0, "M6": 20.0}, ["F3", "F4"], TRIO_B + TRIO_C, [
        T("T01", "M6", "F3", "CR10-60", 575, 3.0e-9, TRIO_C),
        T("T02", "M6", "F3", "CR10-60", 575, 8.6e-9, TRIO_C),
        T("T03", "M3", "F4", "CR6-30", 575, 2.9e-8, TRIO_B),
        T("T04", "M6", "F3", "CR10-60", 600, 3.5e-9, TRIO_C),
        T("T05", "M3", "F4", "CR6-30", 600, 9.9e-9, TRIO_B),
        T("T06", "M3", "F4", "CR6-30", 600, 2.6e-8, TRIO_B),
        T("T07", "M6", "F3", "CR10-60", 625, 3.1e-9, TRIO_C),
        T("T08", "M3", "F4", "CR6-30", 625, 9.4e-9, TRIO_B),
        T("T09", "M3", "F4", "CR6-30", 625, 2.7e-8, TRIO_B),
        T("T10", "M3", "F4", "CR6-30", 625, 6.1e-9, TRIO_B, stop=0.07),
    ], dt.datetime(2023, 11, 6, 10, 15), seed=23)
    # h4: aluminium alloy; a heavier disc on the pan for three tests (log only); a controller excursion late on
    out["h4"] = P("h4", "2618-T61", (188.0, 152.0), {"M2": 10.0, "M5": 10.0}, ["F1", "F5"], TRIO_A + TRIO_B, [
        T("T01", "M2", "F1", "CR8-40", 175, 2.8e-8, TRIO_A),
        T("T02", "M5", "F5", "CR8-40", 175, 9.0e-9, TRIO_B),
        T("T03", "M2", "F1", "CR8-40", 175, 3.3e-9, TRIO_A),
        T("T04", "M5", "F5", "CR8-40", 200, 2.6e-8, TRIO_B, pan_extra=2.00),
        T("T05", "M2", "F1", "CR8-40", 200, 9.2e-9, TRIO_A),
        T("T06", "M5", "F5", "CR8-40", 200, 3.5e-9, TRIO_B, pan_extra=2.00),
        T("T07", "M2", "F1", "CR8-40", 225, 2.5e-8, TRIO_A, excursion=(0.74, 0.80, 9.0)),
        T("T08", "M5", "F5", "CR8-40", 225, 9.5e-9, TRIO_B, pan_extra=2.00),
        T("T09", "M2", "F1", "CR8-40", 225, 3.1e-9, TRIO_A),
    ], dt.datetime(2024, 4, 22, 8, 30), seed=24, notes=[
        ("T08", 2, "Unloading T08: a 2.00 kg disc found under the hanger plate in addition to the masses in the test "
                   "sheet. It has been there since T04 was loaded, so T04, T06 and T08 carried 2.00 kg more than the "
                   "test sheet gives."),
    ])
    # h5: copper on two drawings; thermocouple replaced mid-programme (new certificate), repeated test
    out["h5"] = P("h5", "Cu-OF", (262.0, 39.8), {"M4": 10.0, "M6": 20.0}, ["F6", "F4"], TRIO_A + TRIO_C, [
        T("T01", "M4", "F4", "CR6-30", 250, 2.6e-8, TRIO_A),
        T("T02", "M6", "F6", "CR10-60", 250, 8.7e-9, TRIO_C),
        T("T03", "M4", "F4", "CR6-30", 250, 3.2e-9, TRIO_A),
        T("T04", "M6", "F6", "CR10-60", 275, 2.7e-8, TRIO_C),
        T("T05", "M4", "F4", "CR6-30", 275, 9.1e-9, TRIO_A, swap=(1, "K-149", 0.12, 0.14)),
        T("T06", "M6", "F6", "CR10-60", 275, 3.3e-9, TRIO_C),
        T("T07", "M4", "F4", "CR6-30", 300, 2.5e-8, ("K-104", "K-149", "K-111")),
        T("T08", "M6", "F6", "CR10-60", 300, 9.7e-9, TRIO_C),
        T("T09", "M4", "F4", "CR6-30", 300, 3.4e-9, ("K-104", "K-149", "K-111"), excursion=(0.61, 0.68, -11.0)),
    ], dt.datetime(2024, 9, 9, 9, 15), seed=25, true_ratio={"M6": 19.62}, notes=[
        ("T08", 26, "Lever ratio of M6 checked with a calibrated load cell: 19.62, not the 20.0 in the programme "
                    "sheet. The knife edges had been renewed before T02, so the error applies to T02, T04, T06 and "
                    "T08."),
    ])
    return out
