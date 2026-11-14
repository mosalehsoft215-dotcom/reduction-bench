"""Steady states of the HP and IP sections of a reheat steam turbine during one test point.

Plant model used by the campaign writer. Units: pressure MPa, temperature K, enthalpy kJ/kg,
entropy kJ/(kg K), mass flow kg/s.

  * The main steam flow through the inlet valves (fixed openings, constant throttle pressure)
    follows m = Kv p / sqrt(T).
  * Full-arc units: the packing between the HP and IP inlets is fed from the HP inlet bowl
    (main steam throttled, same enthalpy). Partial-arc units: it is fed from the control-stage
    (first-stage) shell, whose pressure follows the HP blading capacity and whose state follows
    the control-stage efficiency. The leak is a choked labyrinth, m = KL p / sqrt(T), adiabatic.
  * HP blading expands m_ms - m_leak from the HP inlet (or first-stage shell) to the HP exhaust.
  * Reheat spray water (feed-pump interstage) mixes into the cold reheat line.
  * The reheater raises the steam to the hot-reheat temperature with a fixed fractional pressure drop.
  * In the IP inlet bowl the hot reheat steam and the packing leak mix adiabatically; the mixture
    expands through the IP section to the crossover with a fixed section efficiency (defined from the
    IP stop-valve pressure).
  * The IP inlet pressure and the crossover pressure follow flow-capacity laws p = C m sqrt(T).
"""
import math

from iapws import IAPWS97


def h_pT(p, T):
    return IAPWS97(P=p, T=T).h


def s_pT(p, T):
    return IAPWS97(P=p, T=T).s


def s_ph(p, h):
    return IAPWS97(P=p, h=h).s


def T_ph(p, h):
    return IAPWS97(P=p, h=h).T


def h_ps(p, s):
    return IAPWS97(P=p, s=s).h


def expand(p_in, h_in, p_out, eta, p_ref=None):
    """Exit enthalpy of an expansion with isentropic efficiency eta; entropy taken at p_ref (default p_in)."""
    s_in = s_ph(p_ref if p_ref is not None else p_in, h_in)
    hs = h_ps(p_out, s_in)
    return h_in - eta * (h_in - hs)


class Unit:
    def __init__(self, c):
        self.c = c
        n = c["nominal"]
        self.p_ms = n["p_ms"]
        self.Kv = n["m_ms"] * math.sqrt(n["T_ms"]) / self.p_ms
        self._calibrate()

    def _source(self, p_ms, T_ms, h_ms, p_fs):
        c = self.c
        if c["admission"] == "partial-arc":
            h_fs = expand(p_ms, h_ms, p_fs, c["eta_cs"])
            return p_fs, T_ph(p_fs, h_fs), h_fs
        return p_ms * (1 - c["dp_valves"]), T_ms, h_ms

    def _flows_and_states(self, T_ms, T_hrh, m_sp, p_crh, p_hrh, p_x, p_fs, m_leak):
        c = self.c
        p_ms = self.p_ms
        h_ms, s_ms = h_pT(p_ms, T_ms), s_pT(p_ms, T_ms)
        m_ms = self.Kv * p_ms / math.sqrt(T_ms)
        src_p, src_T, src_h = self._source(p_ms, T_ms, h_ms, p_fs)
        if m_leak is None:
            m_leak = self.KL * src_p / math.sqrt(src_T)
        m_hpb = m_ms - m_leak
        h_crh = h_ms - c["eta_hp"] * (h_ms - h_ps(p_crh, s_ms))
        h_w = h_pT(c["spray_p"], c["spray_T"])
        m_hrh = m_hpb + m_sp
        h_crh_mix = (m_hpb * h_crh + m_sp * h_w) / m_hrh
        h_hrh = h_pT(p_hrh, T_hrh)
        m_ip = m_hrh + m_leak
        h_mix = (m_hrh * h_hrh + m_leak * src_h) / m_ip
        h_x = expand(p_hrh, h_mix, p_x, c["eta_ip"])
        return dict(p_ms=p_ms, T_ms=T_ms, h_ms=h_ms, m_ms=m_ms, src_p=src_p, src_T=src_T, leak_h=src_h,
                    m_leak=m_leak, m_hpb=m_hpb, p_crh=p_crh, h_crh=h_crh, T_crh=T_ph(p_crh, h_crh),
                    h_w=h_w, m_sp=m_sp, h_crh_mix=h_crh_mix, T_crh_mix=T_ph(p_crh, h_crh_mix),
                    p_hrh=p_hrh, T_hrh=T_hrh, h_hrh=h_hrh, m_hrh=m_hrh, m_ip=m_ip,
                    h_mix=h_mix, T_mix=T_ph(p_hrh, h_mix), p_x=p_x, h_x=h_x, T_x=T_ph(p_x, h_x),
                    p_fs=p_fs, r=m_leak / m_hrh)

    def _calibrate(self):
        c, n = self.c, self.c["nominal"]
        p_crh, p_x, p_fs = n["p_crh"], n["p_x"], n.get("p_fs", 0.0)
        p_hrh = p_crh * (1 - c["dp_rh"])
        m_leak = c["leak_ratio"] * n["m_ms"] * 0.95
        for _ in range(50):
            st = self._flows_and_states(n["T_ms"], n["T_hrh"], 0.0, p_crh, p_hrh, p_x, p_fs, m_leak)
            new = c["leak_ratio"] * st["m_hrh"]
            if abs(new - m_leak) < 1e-12:
                break
            m_leak = new
        self.KL = m_leak * math.sqrt(st["src_T"]) / st["src_p"]
        self.C_ip = p_hrh / (st["m_ip"] * math.sqrt(st["T_mix"]))
        self.C_lp = p_x / (st["m_ip"] * math.sqrt(st["T_x"]))
        if c["admission"] == "partial-arc":
            self.C_hp = p_fs / (st["m_hpb"] * math.sqrt(st["src_T"]))
        self.nominal_state = st

    def state(self, T_ms, T_hrh, m_sp):
        c, n = self.c, self.c["nominal"]
        st = self.nominal_state
        p_hrh, p_x, p_fs = st["p_hrh"], st["p_x"], st["p_fs"]
        for _ in range(80):
            st = self._flows_and_states(T_ms, T_hrh, m_sp, p_hrh / (1 - c["dp_rh"]), p_hrh, p_x, p_fs, None)
            ph = self.C_ip * st["m_ip"] * math.sqrt(st["T_mix"])
            px = self.C_lp * st["m_ip"] * math.sqrt(st["T_x"])
            pf = self.C_hp * st["m_hpb"] * math.sqrt(st["src_T"]) if c["admission"] == "partial-arc" else p_fs
            d = abs(ph - p_hrh) + abs(px - p_x) + abs(pf - p_fs)
            p_hrh, p_x, p_fs = ph, px, pf
            if d < 1e-11:
                break
        return self._flows_and_states(T_ms, T_hrh, m_sp, p_hrh / (1 - c["dp_rh"]), p_hrh, p_x, p_fs, None)


def apparent_ip(st):
    """IP efficiency computed from hot-reheat conditions as if no steam joined the IP inlet."""
    hs = h_ps(st["p_x"], s_ph(st["p_hrh"], st["h_hrh"]))
    return (st["h_hrh"] - st["h_x"]) / (st["h_hrh"] - hs)
