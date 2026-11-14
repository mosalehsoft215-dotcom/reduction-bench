"""Thermodynamic-method reduction of a boiler-feed-pump test to efficiency and head at the guarantee flow.

python3 bfp_reduce.py <test folder> <answer file>

1. Transport lag: water seen by the discharge thermometers at time t passed the suction thermometers at
   t - tau, tau = (pipe and pump volume between the thermowells) x density / mass flow. Every temperature
   difference pairs a discharge reading with the lagged suction reading.
2. Idle warm-up checks (pump stopped, warm-up flow m_w through it): dT = delta - Q_loss / (m_w cp). Two checks at
   different flows separate the discharge-minus-suction thermometer offset delta from the casing heat loss.
3. Shaft power: switchgear meter through the motor certificate (efficiency linear in load), checked against the
   slip (proportional to shaft output); where the two disagree by more than 1.5 % the slip is taken.
4. Per steady point: E = (p2 - p1)/rho_mean + (c2^2 - c1^2)/2 + g dz; discharge mass flow m from
   P - P_mech - Q_loss = m (h2 - h1 + dKE + g dz) + m_bl (h_bl - h1) with IAPWS-IF97 enthalpies;
   efficiency m E / P; flow and head scaled to rated speed (m n0/n, H (n0/n)^2).
5. Quadratic fits of efficiency and head against flow; values at the guarantee flow.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from iapws import IAPWS97

G = 9.80665
SLIP_TRUST = 0.015


def water(p_bar, **kw):
    if "T" in kw:
        return IAPWS97(P=p_bar / 10.0, T=kw["T"] + 273.15)
    return IAPWS97(P=p_bar / 10.0, h=kw["h"])


def load_test(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "pump_sheet.json").read_text())
    with (folder / "test_record.csv").open() as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], np.array(rows[1:], dtype=float)
    rec = {name: body[:, i] for i, name in enumerate(head)}
    with (folder / "events.csv").open() as fh:
        events = list(csv.DictReader(fh))
    return sheet, rec, events


def windows(events, start, end):
    out, opened = [], None
    for ev in events:
        if ev["event"] == start:
            opened = (float(ev["t_s"]), ev["value"])
        elif ev["event"] == end and opened:
            out.append((opened[0], float(ev["t_s"]), opened[1]))
            opened = None
    return out


def transit_volume(sheet):
    th = sheet["thermometry"]
    s, d = th["suction"]["pipe_to_flange"], th["discharge"]["flange_to_pipe"]
    return (s["length_m"] * math.pi / 4 * s["bore_m"] ** 2 + th["pump_internal_volume_m3"]
            + d["length_m"] * math.pi / 4 * d["bore_m"] ** 2)


def lagged_difference(t, t_suc, t_dis, tau, sel):
    """Mean of discharge minus the suction reading tau seconds earlier, over the selected samples."""
    shifted = np.interp(t - tau, t, t_suc)
    return float(np.mean((t_dis - shifted)[sel])), float(np.mean(shifted[sel]))


def motor_efficiency(sheet, load):
    cert = sheet["motor"]["efficiency_certificate"]
    return float(np.interp(load, [c["load_fraction"] for c in cert], [c["efficiency"] for c in cert]))


def shaft_from_meter(sheet, p_in):
    p_r = sheet["motor"]["rated_output_kw"]
    shaft = p_in * 0.96
    for _ in range(50):
        shaft = p_in * motor_efficiency(sheet, shaft / p_r)
    return shaft


def reduce_test(folder):
    sheet, rec, events = load_test(folder)
    t = rec["t_s"]
    t_suc = 0.5 * (rec["TE-S1A"] + rec["TE-S1B"])
    t_dis = 0.5 * (rec["TE-D1A"] + rec["TE-D1B"])
    vol = transit_volume(sheet)

    # idle warm-up checks: offset and casing heat loss. The downstream thermometer, less the upstream one read
    # tau earlier, is -Q_loss / (m_w cp) plus its own offset relative to the upstream one.
    reverse = sheet["warm_up_line"]["direction"] == "discharge to suction"
    up, down = (t_dis, t_suc) if reverse else (t_suc, t_dis)
    checks = []
    for a, b, _ in windows(events, "warm-up check start", "warm-up check end"):
        sel = (t >= a + 0.4 * (b - a)) & (t < b)
        p1 = float(np.mean(rec["PT-S"][sel]))
        props = water(p1, T=float(np.mean(up[sel])))
        flow = rec["FT-WU"]
        # plug flow: the delay is the time over which the warm-up flow has swept the transit volume
        swept = np.concatenate([[0.0], np.cumsum(0.5 * (flow[1:] + flow[:-1]) / props.rho * np.diff(t))])
        tau = t - np.interp(swept - vol, swept, t)
        d, _ = lagged_difference(t, up, down, tau, sel)
        checks.append((float(np.mean(flow[sel])), d, props.cp))
    (m_a, d_a, cp_a), (m_b, d_b, cp_b) = checks[0], checks[-1]
    # d = offset(down - up) - Q x, x = 1 / (m cp): two equations
    x_a, x_b = 1.0 / (m_a * cp_a), 1.0 / (m_b * cp_b)
    q_loss = (d_b - d_a) / (x_a - x_b)
    offset = d_a + q_loss * x_a
    delta = -offset if reverse else offset  # discharge minus suction

    n0 = sheet["rated_speed_rpm"]
    p_r = sheet["motor"]["rated_output_kw"]
    s_r = 1.0 - sheet["motor"]["rated_speed_rpm_at_50hz"] / 3000.0
    a_s = math.pi / 4 * sheet["suction_tapping"]["bore_m"] ** 2
    a_d = math.pi / 4 * sheet["discharge_tapping"]["bore_m"] ** 2
    dz = sheet["discharge_tapping"]["elevation_m"] - sheet["suction_tapping"]["elevation_m"]
    pm_r = sheet["bearing_and_seal_losses"]["kw_at_rated_speed"]

    points = []
    for a, b, _ in windows(events, "steady start", "steady end"):
        sel = (t >= a) & (t < b)
        mean = {k: float(np.mean(v[sel])) for k, v in rec.items()}
        n, f = mean["SE-P"], mean["F-HZ"]
        slip = 1.0 - n / (60.0 * f)
        points.append(dict(sel=sel, mean=mean, n=n, meter=shaft_from_meter(sheet, mean["JM-1"]),
                           by_slip=p_r * slip / s_r))
    ratio = float(np.median([pt["by_slip"] / pt["meter"] for pt in points]))
    use_slip = abs(ratio - 1.0) > SLIP_TRUST

    flows, heads, effs = [], [], []
    for pt in points:
        mean, sel, n = pt["mean"], pt["sel"], pt["n"]
        shaft = pt["by_slip"] if use_slip else pt["meter"]
        p1, p2 = mean["PT-S"], mean["PT-D"]
        p_mech = pm_r * (n / n0) ** 2
        h_bl = water(mean["PT-BL"], T=mean["TE-BL"]).h
        m = 0.9 * shaft / 20.0
        for _ in range(30):
            rho_guess = water(p1, T=float(np.mean(t_suc[sel]))).rho
            tau = vol * rho_guess / m
            d, t1 = lagged_difference(t, t_suc, t_dis, tau, sel)
            t2 = t1 + d - delta
            s1, s2 = water(p1, T=t1), water(p2, T=t2)
            rbar = 0.5 * (s1.rho + s2.rho)
            c1, c2 = m / (s1.rho * a_s), m / (s2.rho * a_d)
            extra = ((c2 * c2 - c1 * c1) / 2 + G * dz) / 1000.0
            m_new = (shaft - p_mech - q_loss - mean["FT-BL"] * (h_bl - s1.h)) / (s2.h - s1.h + extra)
            if abs(m_new - m) < 1e-7:
                m = m_new
                break
            m = m_new
        e = (p2 - p1) * 1e5 / rbar + (c2 * c2 - c1 * c1) / 2 + G * dz
        flows.append(m * n0 / n)
        heads.append(e / G * (n0 / n) ** 2)
        effs.append(m * e / (shaft * 1000.0))
    m_g = sheet["guarantee_flow_kg_s"]
    eta = float(np.polyval(np.polyfit(flows, effs, 2), m_g))
    head = float(np.polyval(np.polyfit(flows, heads, 2), m_g))
    notes = dict(delta=delta, q_loss=q_loss, slip_over_meter=ratio, used="slip" if use_slip else "meter",
                 flows=[round(x, 3) for x in flows], effs=[round(x, 5) for x in effs])
    return {"efficiency_pct_at_guarantee": round(100.0 * eta, 4), "head_m_at_guarantee": round(head, 3)}, notes


def main():
    answer, notes = reduce_test(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(answer, indent=1) + "\n")
    print(json.dumps(notes))


if __name__ == "__main__":
    main()
