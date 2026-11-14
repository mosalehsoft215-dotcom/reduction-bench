"""Water-cooled stator of generator G-20: winding layout, bar losses, bar hydraulics and steady temperatures.

All 84 bars (42 slots x top/bottom) are in parallel between the inlet and outlet water headers. A bar's flow
follows dp = K q^2. Its losses are I_bar^2 (R_dc(T) + R_eddy(T)), with R_dc rising and R_eddy falling with copper
temperature, and all of them leave in the bar's own water. Copper sits above the mean water temperature by
loss / G_film. A slot RTD reads the mean copper temperature of the two bars in its slot.
"""
import math

SLOTS, PATHS, PITCH = 42, 2, 17
RATED_KV, RATED_MVA = 20.0, 600.0
CP_WATER, ALPHA_CU = 4180.0, 0.00393
BELTS = ["A", "C", "B", "A", "C", "B"]
K_TOP, K_BOTTOM = 1.662e6, 1.50e6
RDC20 = 1.97e-4
EDDY_SHARE = {"top": 0.45, "bottom": 0.15}
G_FILM = 3000.0
HEAT_CAP = {"top": 2.9e5, "bottom": 3.0e5}
SINGLE_FEED_GRADIENT = 0.18


def bars():
    out = []
    for layer in ("top", "bottom"):
        for s in range(1, SLOTS + 1):
            top_slot = s if layer == "top" else (s - 1 - PITCH) % SLOTS + 1
            out.append({
                "bar": f"{'T' if layer == 'top' else 'B'}{s:02d}",
                "slot": s,
                "layer": layer,
                "phase": BELTS[(top_slot - 1) // 7],
                "junction_box": f"JB-{(s - 1) // 7 + 1}",
                "k_pa_per_kg2_s2": K_TOP if layer == "top" else K_BOTTOM,
                "r_dc_20c_ohm": RDC20,
                "r_eddy_20c_ohm": RDC20 * EDDY_SHARE[layer],
                "g_film_w_per_k": G_FILM,
                "heat_capacity_j_per_k": HEAT_CAP[layer],
                "header_position": round((s - 0.5) / SLOTS, 4),
            })
    return out


def bar_current(phase_amps, phase):
    return phase_amps[phase] / PATHS


def loss_w(i_bar, bar, t_cu):
    x = 1.0 + ALPHA_CU * (t_cu - 20.0)
    return i_bar * i_bar * (bar["r_dc_20c_ohm"] * x + bar["r_eddy_20c_ohm"] / x)


def steady_bar(bar, i_bar, t_in, q):
    """Outlet water and copper temperature of one bar at flow q (kg/s)."""
    t_cu = t_in + 10.0
    for _ in range(60):
        loss = loss_w(i_bar, bar, t_cu)
        t_out = t_in + loss / (q * CP_WATER)
        t_new = 0.5 * (t_in + t_out) + loss / bar["g_film_w_per_k"]
        if abs(t_new - t_cu) < 1e-10:
            break
        t_cu = t_new
    return t_out, t_cu, loss


def bar_flow(bar, dp_pa, fraction=1.0):
    return fraction * math.sqrt(dp_pa / bar["k_pa_per_kg2_s2"])


def local_dp(bar, dp_header, single_feed):
    """dp across a bar: the header dp with both feeds open; with the exciter-end feed only it falls along the header."""
    return dp_header * (1.0 - SINGLE_FEED_GRADIENT * bar["header_position"]) if single_feed else dp_header


def time_constant_s(bar, q):
    return bar["heat_capacity_j_per_k"] / (q * CP_WATER)
