# lunar_tes.py -- sintered-regolith TES as N brick segments in a row, with two gas loops in
# counterflow (charging A -> B with sun, discharging B -> A always).
# Loaded by Step 3 of the notebook, after lunar_regolith.py (keep this file next to the notebook).
# Not a stand-alone script: it uses np, CP, say() (Step 1), the inputs (Step 2) and
# cp_regolith(), e_regolith(), k_regolith() (lunar_regolith.py).
# =====================================================================
#   Segments 0 .. N-1. Segment 0 = END A (hot end), segment N-1 = END B (cold end).
#   Each segment: brick mass m_TES / N at one temperature T_i.
#   Charging gas (from the receiver, T_rec_set) passes segments 0 -> N-1.
#   Discharging gas (back from the PHX) passes segments N-1 -> 0 and leaves END A to the PHX.
#   Perfect heat exchange in each segment: the gas leaves segment i at T_i, so
#       heat into segment i = m_gas * (h(T_gas entering i) - h(T_i))      (> 0 charges it)
#   With several segments in a row this gives a hot end, a cold end and a temperature front
#   between them; more segments = sharper front.
#   Heat loss: the store is a cube of volume m_TES / rho_TES under L_ins of loose regolith;
#   each segment loses 1/N of the cube surface: k_reg(T_avg) * A / L_ins * (T_i - T_ground).
#   Pressure drop per loop: N_segments * dp_seg + dp_other_frac * P_loop.
#   advance() moves the segments over dt hours; it splits dt into sub-steps so no segment can
#   take more than ~30 % of the heat that would bring it to the gas temperature in one sub-step.
# =====================================================================

# --- sCO2 at P_loop (gas side) -------------------------------------------------------------
_Tg = np.linspace(300.0, T_rec_set + 60.0, 800)
_hg = np.array([CP.PropsSI('H', 'T', T, 'P', P_loop, FLUID) for T in _Tg])
_rhog = np.array([CP.PropsSI('D', 'T', T, 'P', P_loop, FLUID) for T in _Tg])
_cpg = np.array([CP.PropsSI('C', 'T', T, 'P', P_loop, FLUID) for T in _Tg])

def gas_h(T):   return np.interp(T, _Tg, _hg)          # J/kg
def gas_T(h):   return np.interp(h, _hg, _Tg)          # K
def gas_rho(T): return np.interp(T, _Tg, _rhog)        # kg/m3

# --- brick energy <-> temperature (segments at end B can sit below T_TES_min) ---------------
_Tb = np.linspace(360.0, T_rec_set + 20.0, 3000)
_eb = e_regolith(_Tb)                                   # J/kg, relative to T_TES_min

def brick_T(e): return np.interp(e, _eb, _Tb)

# --- store --------------------------------------------------------------------------------
def tes_geometry(m_TES):
    """Segment mass and insulation area for a TES (brick) mass m_TES, kg."""
    side = (m_TES / rho_TES) ** (1 / 3)
    return dict(m_TES=m_TES, m_seg=m_TES / N_segments, side=side,
                A_loss=6 * side**2 / N_segments)

def _pass(T_s, order, m_gas, T_in):
    """Gas through the segments in `order`: heat to each segment (W) and gas outlet (K)."""
    Q = np.zeros(N_segments)
    if m_gas <= 0:
        return Q, T_in
    h_in = gas_h(T_in)
    for i in order:
        h_i = gas_h(T_s[i])
        Q[i] = m_gas * (h_in - h_i)
        h_in = h_i                                       # gas leaves at the segment temperature
    return Q, T_s[order[-1]]

ORDER_CHARGE = np.arange(N_segments)                    # A -> B
ORDER_DISCHARGE = ORDER_CHARGE[::-1]                    # B -> A

def charge_pass(T_s, m_gas, T_in=None):
    """Charging gas (receiver outlet, default T_rec_set) through the segments A -> B."""
    return _pass(T_s, ORDER_CHARGE, m_gas, T_rec_set if T_in is None else T_in)

def discharge_pass(T_s, m_gas, T_in):
    """Discharging gas (back from the PHX at T_in) through the segments B -> A."""
    return _pass(T_s, ORDER_DISCHARGE, m_gas, T_in)

def loop_dp():
    """Loop pressure drop, Pa: bricks (dp_seg per segment) + receiver or PHX + piping."""
    return N_segments * dp_seg + dp_other_frac * P_loop

def blower_power(m_gas, T_inlet):
    """Blower electric power, kWe (blower on the cold leg, inlet temperature T_inlet)."""
    return m_gas * loop_dp() / (gas_rho(T_inlet) * eta_blower) / 1e3 if m_gas > 0 else 0.0

def heat_loss(T_s, geo):
    """Heat lost through the insulation from each segment, W."""
    return k_regolith(0.5 * (T_s + T_ground)) * geo['A_loss'] / L_ins * (T_s - T_ground)

def stored_energy(T_s, geo):
    """Energy stored above T_TES_min, kWh."""
    return float(np.sum(e_regolith(T_s)) * geo['m_seg'] / 3.6e6)

def advance(T_s, dt_h, geo, m_ch=0.0, T_ch_in=None, m_dis=0.0, T_dis_in=None):
    """Move the segments over dt_h hours with fixed flows (kg/s) and inlet temperatures.
    Returns new segment temperatures and step averages: heat into the bricks from the charging
    gas, heat taken by the discharging gas, heat loss (kW), and the gas outlet temperatures
    (end B for charging, end A for discharging, K)."""
    T_s = np.array(T_s, dtype=float)
    cap = geo['m_seg'] * cp_regolith(T_s.min())        # J/K per segment
    flow_cap = (m_ch * np.interp(T_rec_set, _Tg, _cpg) if m_ch > 0 else 0.0) + \
               (m_dis * np.interp(T_s.max(), _Tg, _cpg) if m_dis > 0 else 0.0)
    n_sub = max(1, int(np.ceil(flow_cap * dt_h * 3600 / (0.3 * cap))))
    dts = dt_h * 3600 / n_sub
    acc = dict(Q_ch=0.0, Q_dis=0.0, Q_loss=0.0, T_ch_out=0.0, T_dis_out=0.0)
    for _ in range(n_sub):
        Qc, Tc = charge_pass(T_s, m_ch, T_ch_in) if m_ch > 0 else (np.zeros(N_segments), np.nan)
        Qd, Td = discharge_pass(T_s, m_dis, T_dis_in) if m_dis > 0 else (np.zeros(N_segments), np.nan)
        Ql = heat_loss(T_s, geo)
        T_s = brick_T(e_regolith(T_s) + (Qc + Qd - Ql) * dts / geo['m_seg'])
        acc['Q_ch'] += Qc.sum(); acc['Q_dis'] -= Qd.sum(); acc['Q_loss'] += Ql.sum()
        acc['T_ch_out'] += Tc; acc['T_dis_out'] += Td
    out = {k: v / n_sub for k, v in acc.items()}
    for k in ('Q_ch', 'Q_dis', 'Q_loss'):
        out[k] /= 1e3                                    # kW
    out['n_sub'] = n_sub
    return T_s, out

_g = tes_geometry(1.0e6)
say(f"TES: {N_segments} brick segments; per 1000 t a cube of {_g['side']:.1f} m side, "
    f"{_g['m_seg']/1e3:.0f} t per segment | blower dp {loop_dp()/1e3:.0f} kPa per loop")