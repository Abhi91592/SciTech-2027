# lunar_plant.py -- plant model of the lunar regolith-TES sCO2 system (layout of paper Fig. 2)
# Loaded by Step 6 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses the inputs (Step 2), the regolith properties (Step 3),
# the cycle maps (Steps 4-5), the solar functions (Step 8 is loaded later; T_surface and
# T_sink_t are only called once the hourly runs start) and np, CP, say() from Step 1.
# ---------------------------------------------------------------------
# Plant performance at one TES temperature (and time), for cycle ct -- layout of Fig. 2
#   t     = hours since sunrise (sets day/night loads and the radiator sink, Step 8); None = night
#   Q_rec = solar heat the receiver puts into the TES gas loop this hour (kW, Step 9); 0 = night
#
# TES gas loop (sCO2 at P_loop): receiver -> TES -> PHX -> gas tank -> blower -> receiver
#   Gas leaving the TES (= PHX gas inlet): T_TES + dT_TES_HX while charging, T_TES - dT_TES_HX
#   otherwise. PHX (counterflow, gas -> cycle sCO2): TIT = gas inlet - dT_PHX. The gas flow
#   is the smallest that keeps every point of the PHX >= dT_PHX apart (r_gas = m_gas/m_cycle,
#   pre-computed on the cycle map). By day it is also >= the flow that carries the heat
#   stored in the TES with a rise of dT_charge: m_gas >= (Q_rec - Q_in) / (h(T_ex + dT_charge) - h(T_ex)).
#   Gas leaving the PHX: h_out = h_in - m_cycle * q_in / m_gas; it returns through the tank
#   and the blower to the receiver (by day) or straight to the TES (night, receiver bypassed).
#   Receiver: gas from h_out to h_out + Q_rec / m_gas; T_rec = mean gas temperature + dT_wall.
#   By day the gas flow is capped at the loop's design flow; solar heat above what that flow
#   can carry (Q_rec_max) is curtailed (heliostats defocused, Step 9).
#   Blower: W = m_gas * dp / (rho(blower inlet) * eta_blower). The loop pressure drop grows
#   with flow like any duct: dp = dp_loop * P_loop * (m_gas^2/rho) / (m_gas^2/rho)_design, where
#   the design point is the largest m_gas^2/rho of normal operation (found in Step 10;
#   until then dp = dp_loop * P_loop).
# Power cycle: tables from Steps 4-5 (state 1 = compressor inlet, after the precooler).
# Coolant loop (water at P_cool_loop): precooler -> habitat heat exchanger -> radiator -> pump
#   Precooler (counterflow, cycle sCO2 -> water): water enters at T1 - dT_precool; its flow is
#   the smallest that keeps every point >= dT_precool apart and the water outlet <= T_cool_max
#   (r_w = m_water/m_cycle, pre-computed on the cycle map). The habitat takes the heat it
#   needs from the hot water first; the rest goes to the radiator panels, whose temperature
#   falls with the water. Area per kg/s of water:
#       A/m_w = integral of dh / (rad_sides*emissivity*sigma*(T(h)^4 - T_sink^4))
#   The radiator is fixed hardware (A_RAD, sized at night). Each hour the plant uses the
#   COLDEST compressor inlet T1 >= T_comp_in the radiator allows. If even T_comp_in_max is
#   not enough, the cycle runs at REDUCED POWER (flow cut until the waste heat fits).
#   Pump: W = m_water * dp_cool_loop / (rho_water * eta_pump).
# Parasitic loads: W_aux = W_aux_fixed + heliostat drives (W_extra, day) + blower + pump;
#   W_net = (W_elec_target + W_aux) / eta_gen, solved by a short fixed-point loop
#   (blower and pump grow with the cycle flow).
# Habitat (Step 2): each hour
#   wall loss   Q_wall = eps_MLI*sigma*A_hab*(T_hab^4 - T_outside^4)   (< 0 = gain by day)
#   inside heat Q_int  = f_elec_inside*(W_elec_target - W_hp) + crew metabolic
#   Q_int < Q_wall -> HEATING need, supplied from the coolant-loop water
#   Q_int > Q_wall -> COOLING need, rejected by the habitat's own radiator at
#                     T_hab_rad, or with a heat pump (COP) to T_hp_cond when the sky
#                     is too warm. The heat pump runs on the habitat's 100 kWe.
# ---------------------------------------------------------------------
# --- Property tables (linear interpolation; h(T) and T(h) are exact inverses on a table) ---
def _table(fluid, P, T_lo, T_hi, n):
    T = np.linspace(T_lo, T_hi, n)
    return T, np.array([CP.PropsSI('H', 'T', x, 'P', P, fluid) for x in T]), \
        np.array([CP.PropsSI('D', 'T', x, 'P', P, fluid) for x in T])

_T_pc, _h_pc, _ = _table(FLUID, P_low, T_comp_in - 5, TIT_top + 5, 1500)        # cycle sCO2, precooler
_T_g, _h_g, _rho_g = _table(FLUID, P_loop, 260.0, T_charge_out_max + 10, 1500)  # TES gas loop
_T_w, _h_w, _rho_w = _table(COOLANT, P_cool_loop, T_comp_in - dT_precool - 5, T_cool_max + 5, 600)
_rho_w_pump = float(np.interp(T_comp_in - dT_precool, _T_w, _rho_w))   # pump on the cold water
_T_phx = np.linspace(T_comp_in, TIT_top + 0.5, 1500)                    # cycle sCO2 in the PHX
_h_phx = {ct: np.array([CP.PropsSI('H', 'T', T, 'P', P_phx_out[ct], FLUID) for T in _T_phx])
          for ct in cycles}
_f = np.linspace(0.0, 1.0, 61)                                          # points along an exchanger

_rad_G_cache = {}                       # G(h) for each sink temperature already used

def _rad_G(T_sink):
    """G(h) on the water table = radiator area per kg/s of water to cool from the top of the
    table down to h, m2/(kg/s). Area from h_a (hot) down to h_b (cold) = G(h_b) - G(h_a).
    Calculated once per sink temperature and reused (same numbers)."""
    key = float(T_sink)
    if key not in _rad_G_cache:
        if len(_rad_G_cache) > 5000:
            _rad_G_cache.clear()
        _rad_G_cache[key] = _rad_G_calc(key)
    return _rad_G_cache[key]

def _rad_G_calc(T_sink):
    with np.errstate(divide='ignore'):
        f = np.where(_T_w > T_sink + 1e-6,
                     1.0 / (rad_sides * emissivity * sigma * (_T_w**4 - T_sink**4)), np.inf)
    seg = 0.5 * (f[1:] + f[:-1]) * np.diff(_h_w)
    return np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])

def _bisect_ratio(ok_fn, n, lo=1e-3, hi=1e3, it=50):
    """Smallest ratio (vector of n) with ok_fn(ratio) True; ok_fn is monotonic (log bisection)."""
    lo, hi = np.full(n, lo), np.full(n, hi)
    for _ in range(it):
        mid = np.sqrt(lo * hi)
        ok = ok_fn(mid)
        hi, lo = np.where(ok, mid, hi), np.where(ok, lo, mid)
    return hi

def phx_ratio(TIT, h_c_in, ct):
    """Smallest m_gas/m_cycle in the PHX, gas inlet TIT + dT_PHX (vector over h_c_in)."""
    h_c_in = np.atleast_1d(np.asarray(h_c_in, dtype=float))
    out = np.full(h_c_in.shape, np.nan)
    v = np.isfinite(h_c_in)
    if not v.any():
        return out
    h_c_out = np.interp(TIT, _T_phx, _h_phx[ct])
    hc = h_c_in[v, None] + _f[None, :] * (h_c_out - h_c_in[v, None])     # cycle side, cold -> hot
    Tc = np.interp(hc, _h_phx[ct], _T_phx)
    hg_in = np.interp(TIT + dT_PHX, _T_g, _h_g)
    def ok(r):
        Tg = np.interp(hg_in - (h_c_out - hc) / r[:, None], _h_g, _T_g)   # gas at the same place
        return np.min(Tg - Tc, axis=1) >= dT_PHX - 1e-4
    out[v] = _bisect_ratio(ok, int(v.sum()))
    return out

def precooler_ratio(h_cool_in, h1):
    """Smallest m_water/m_cycle in the precooler (vector): every point >= dT_precool apart,
    water outlet <= T_cool_max. Water enters at T1 - dT_precool."""
    h_cool_in, h1 = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (h_cool_in, h1))
    out = np.full(h1.shape, np.nan)
    v = np.isfinite(h_cool_in) & np.isfinite(h1)
    if not v.any():
        return out
    a, b = h1[v, None], h_cool_in[v, None]
    T1 = np.interp(h1[v], _h_pc, _T_pc)
    hw_c = np.interp(T1 - dT_precool, _T_w, _h_w)[:, None]
    hw_max = np.interp(T_cool_max, _T_w, _h_w)
    Tco2 = np.interp(a + _f[None, :] * (b - a), _h_pc, _T_pc)                 # cold end -> hot end
    def ok(r):
        hw = hw_c + _f[None, :] * (b - a) / r[:, None]
        Tw = np.interp(hw, _h_w, _T_w)
        return (np.min(Tco2 - Tw, axis=1) >= dT_precool - 1e-4) & (hw[:, -1] <= hw_max + 1e-6)
    out[v] = _bisect_ratio(ok, int(v.sum()))
    return out

# --- Add the loop flow ratios to the cycle maps (Steps 4-5), so each hour is a table lookup ---
for ct in cycles:
    G = CYC[ct]['G']
    G['r_gas'] = np.column_stack([phx_ratio(TIT_grid[j], G['h_phx_in'][:, j], ct)
                                  for j in range(len(TIT_grid))])
    G['r_w'] = precooler_ratio(G['h_cool_in'], G['h1'])
KEYS = KEYS + ['r_gas', 'r_w']
def _ratio_range(ct, k):
    G = CYC[ct]['G']
    with np.errstate(invalid='ignore'):
        ok = (G['w_net'] > 0) & (G['q_in'] > 0) & (G['eta'] >= eta_cycle_min)   # working points
    return f"{ct} {np.nanmin(G[k][ok]):.2f}-{np.nanmax(G[k][ok]):.2f}"
say("Loop flow per kg/s of cycle sCO2: TES gas " + ", ".join(_ratio_range(ct, 'r_gas') for ct in cycles)
    + "; water " + ", ".join(_ratio_range(ct, 'r_w') for ct in cycles))

def habitat_balance(t=None):
    """Habitat heat balance at time t (None = night). All in kW."""
    T_out = T_surf_night if t is None else float(T_surface(t))
    T_sky = T_sink_night if t is None else float(T_sink_t(t))
    Q_wall = eps_MLI * sigma * A_hab * (T_hab**4 - T_out**4) / 1e3   # + lost, - gained
    Q_met = n_crew * q_metab
    hp_on = T_sky > T_hab_rad - dT_hab_rad                         # sky too warm: heat pump
    COP = f_carnot_hp * T_hp_evap / (T_hp_cond - T_hp_evap)
    if hp_on:   # Q_net = f*(W - Q_net/COP) + Q_met - Q_wall, solved for Q_net
        Q_net = (f_elec_inside * W_elec_target + Q_met - Q_wall) / (1 + f_elec_inside / COP)
    else:
        Q_net = f_elec_inside * W_elec_target + Q_met - Q_wall
    Q_cool = max(Q_net, 0.0)
    Q_heat = max(-Q_net, 0.0)
    W_hp = Q_cool / COP if hp_on else 0.0
    T_r = T_hp_cond if hp_on else T_hab_rad
    q_r = rad_sides * emissivity * sigma * (T_r**4 - T_sky**4) / 1e3  # kW per m2 of panel
    Q_hrad = Q_cool + W_hp                                          # heat to habitat radiator
    return dict(Q_wall=Q_wall, Q_int=f_elec_inside * (W_elec_target - W_hp) + Q_met,
                Q_heat=Q_heat, Q_cool=Q_cool, W_hp=W_hp, hp_on=hp_on, Q_hrad=Q_hrad,
                A_hrad=Q_hrad / q_r if q_r > 0 else np.inf, T_out=T_out)

def _loop_at_flow(c, m_c, L):
    """Both loops for cycle values c (scalars or arrays over T1) at cycle flow m_c (kg/s)."""
    m_c = np.asarray(m_c, dtype=float)
    Q_in = m_c * c['q_in'] / 1e3                                         # kW to the cycle
    m_ch = np.maximum(L['Q_rec'] - Q_in, 0.0) * 1e3 / L['dh_ch']          # kg/s to carry the stored heat
    # By day the gas flow may not exceed the loop's design flow (m_gas^2/rho = design): above
    # it the heliostats are defocused. m_cap solves m^2 / rho(blower inlet at that flow) = ref.
    # (Only needed while the sun charges the loop: at night Q_rec = 0, so no gas flow is
    #  set by charging and the cap cannot act -- it is skipped, results are unchanged.)
    ref = GAS_DP_REF[L['ct']]
    with np.errstate(invalid='ignore', divide='ignore'):
        if ref and L['Q_rec'] > 0:
            m_cap = np.sqrt(ref * np.interp(L['hg_in'], _h_g, _rho_g)) * np.ones_like(Q_in)
            for _ in range(6):
                m_cap = np.sqrt(ref * np.interp(L['hg_in'] - Q_in * 1e3 / m_cap, _h_g, _rho_g))
        else:
            m_cap = np.full_like(Q_in, np.inf)
    Q_rec_max = Q_in + m_cap * L['dh_ch'] / 1e3                         # kW the loop can take in
    m_ch = np.minimum(m_ch, m_cap)
    m_g = np.maximum(c['r_gas'] * m_c, m_ch)
    with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
        hg_out = np.where(m_g > 0, L['hg_in'] - m_c * c['q_in'] / np.where(m_g > 0, m_g, 1.0), L['hg_in'])
        rho_g = np.interp(hg_out, _h_g, _rho_g)                            # blower inlet (coldest gas)
        dp = dp_loop * P_loop * (m_g**2 / rho_g / ref if ref else 1.0)   # Pa, loop pressure drop
        W_bl = m_g * dp / (rho_g * eta_blower) / 1e3                      # kWe
        m_w = c['r_w'] * m_c
        W_pump = m_w * dp_cool_loop / (_rho_w_pump * eta_pump) / 1e3     # kWe
        Q_r = np.minimum(L['Q_rec'], Q_rec_max)                          # heat the receiver passes on
        h_ro = np.where(m_g > 0, hg_out + Q_r * 1e3 / np.where(m_g > 0, m_g, 1.0), hg_out)
    T_g_out = np.interp(hg_out, _h_g, _T_g)
    T_ro = np.interp(h_ro, _h_g, _T_g)
    T_rec = np.where(L['Q_rec'] > 0, 0.5 * (T_g_out + T_ro) + dT_wall, np.nan)
    return dict(m_g=m_g, T_g_out=T_g_out, T_ro=T_ro, T_rec=T_rec, W_bl=W_bl, m_w=m_w, W_pump=W_pump,
                dp=dp, rho_g=rho_g, Q_rec_max=Q_rec_max)

def _loop(c, L):
    """Shaft power for the habitat + all parasitic loads, with the loops at that flow."""
    W_net = np.full_like(np.asarray(c['w_net'], dtype=float), L['W_fix'] / eta_gen)
    with np.errstate(invalid='ignore', over='ignore'):
        for _ in range(40):
            lp = _loop_at_flow(c, W_net * 1e3 / c['w_net'], L)
            W_new = (L['W_fix'] + lp['W_bl'] + lp['W_pump']) / eta_gen
            diff = np.abs(W_new - W_net)
            W_net = W_new
            if not np.any(diff > 1e-7):
                break
        # no solution where the blower + pump grow faster than the power they cost
        # (very high flow near the TES floor): those points cannot give full power
        W_net = np.where(np.isfinite(W_net) & (diff <= 1e-6) & (W_net < 10 * L['W_fix'] / eta_gen),
                         W_net, np.nan)
        lp = _loop_at_flow(c, W_net * 1e3 / c['w_net'], L)
    lp['W_net'] = W_net
    return lp

def _area_needed(c, m_w, T_sink, Q_heat=0.0):
    """Radiator area needed (m2) for cycle values c and water flow m_w (arrays or scalars)."""
    T1 = np.interp(c['h1'], _h_pc, _T_pc)
    hw_c = np.interp(T1 - dT_precool, _T_w, _h_w)
    dh_w = c['q_out'] / c['r_w']                                         # water enthalpy rise
    with np.errstate(invalid='ignore', divide='ignore'):
        f_h = np.clip(Q_heat * 1e3 / (m_w * dh_w), 0.0, 1.0)              # share to the habitat
        h_rad_in = hw_c + (1 - f_h) * dh_w                                # after the habitat HX
        G = _rad_G(T_sink)
        return m_w * (np.interp(hw_c, _h_w, G) - np.interp(h_rad_in, _h_w, G))

A_RAD = {ct: np.inf for ct in cycles}     # sized below, once the plant model exists

GAS_DP_REF = {ct: None for ct in cycles}   # design m_gas^2/rho of the gas loop (Step 10)

def _loop_context(T_TES, W_extra, Q_rec, ct):
    """TES gas exit (= PHX gas inlet) and the fixed parts of the loop calculation."""
    charging = Q_rec > 0
    T_ex = T_TES + dT_TES_HX if charging else T_TES - dT_TES_HX
    hg_in = float(np.interp(T_ex, _T_g, _h_g))
    return dict(T_ex=T_ex, hg_in=hg_in, Q_rec=float(Q_rec), W_extra=float(W_extra), ct=ct,
                dh_ch=float(np.interp(T_ex + dT_charge, _T_g, _h_g)) - hg_in,
                W_fix=W_elec_target + W_aux_fixed + W_extra)

def plant_off(T_TES, t=None, ct='simple', L=None):
    """Plant state when the cycle is off: no heat to the cycle, no electricity. By day the
    blower still moves the gas that charges the TES."""
    L = L or _loop_context(T_TES, 0.0, 0.0, ct)
    m_g = L['Q_rec'] * 1e3 / L['dh_ch'] if L['Q_rec'] > 0 else 0.0
    rho = float(np.interp(L['hg_in'], _h_g, _rho_g))
    ref = GAS_DP_REF[ct]
    dp = dp_loop * P_loop * (m_g**2 / rho / ref if ref else 1.0)
    W_bl = m_g * dp / (rho * eta_blower) / 1e3
    return dict(T_TES=T_TES, TIT=np.nan, dT_PHX_t=np.nan, h_phx_in=np.nan, too_cold=False,
                eta=np.nan, m_dot=0.0,
                Q_in=0.0, Q_out=0.0, W_c=0.0, W_t=0.0, W_net=0.0, Q_recup=0.0, W_gross=0.0,
                W_aux=W_aux_fixed + L['W_extra'] + W_bl, W_elec=0.0, Q_hab=0.0,
                Q_rad=0.0, A_rad=A_RAD[ct], A_need=np.inf, rad_use=np.nan,
                T_sink=T_sink_night if t is None else float(T_sink_t(t)),
                T_cooler_in=np.nan, T_comp_in=np.nan, running=False, part_load=False, load=0.0,
                m_gas=m_g, W_blower=W_bl, W_pump=0.0, m_water=0.0, T_gas_hot=L['T_ex'],
                dp_loop_t=dp / P_loop, K_gas=m_g**2 / rho, Q_rec_max=np.inf,
                T_gas_cold=np.nan, T_rec=np.nan, T_rec_out=np.nan, T_w_hot=np.nan, T_w_cold=np.nan,
                **_hab_fields(habitat_balance(t), 0.0))

def _hab_fields(hb, Q_hab):
    """Habitat results to record, given the heat the cycle actually delivered."""
    return dict(Q_hab_need=hb['Q_heat'], Q_hab_unmet=max(hb['Q_heat'] - Q_hab, 0.0),
                Q_hab_cool=hb['Q_cool'], W_hp=hb['W_hp'], hp_on=hb['hp_on'],
                Q_hab_wall=hb['Q_wall'], Q_hab_rad=hb['Q_hrad'], A_hab_rad=hb['A_hrad'])

def _plant_state(T_TES, TIT, c, T1, W_net, A_need, T_sink, ct, load, hb, L):
    """All plant flows for cycle point c (scalars) at shaft power W_net."""
    m_dot = W_net * 1e3 / c['w_net']                      # kg/s cycle sCO2
    lp = _loop_at_flow(c, m_dot, L)
    Q_in = m_dot * c['q_in'] / 1e3                        # kW from the TES gas (PHX)
    Q_out = m_dot * c['q_out'] / 1e3                      # kW waste heat (precooler)
    Q_hab = min(hb['Q_heat'], Q_out)                      # heat the habitat needs (Eq. 26)
    Q_rad = Q_out - Q_hab                                 # rest to the radiator (Eq. 27)
    W_aux = W_aux_fixed + L['W_extra'] + float(lp['W_bl']) + float(lp['W_pump'])   # kWe
    W_gross = W_net * eta_gen                             # kWe at the generator
    T_wc = T1 - dT_precool
    T_wh = float(np.interp(np.interp(T_wc, _T_w, _h_w) + c['q_out'] / c['r_w'], _h_w, _T_w))
    return dict(T_TES=T_TES, TIT=TIT, dT_PHX_t=T_TES - TIT, h_phx_in=c['h_phx_in'], too_cold=False,
                eta=c['eta'], m_dot=m_dot, Q_in=Q_in, Q_out=Q_out,
                W_c=m_dot * c['w_c'] / 1e3, W_t=m_dot * c['w_t'] / 1e3, W_net=W_net,
                Q_recup=m_dot * c['q_rec'] / 1e3, W_gross=W_gross, W_aux=W_aux,
                W_elec=W_gross - W_aux, Q_hab=Q_hab, Q_rad=Q_rad,
                A_rad=A_RAD[ct], A_need=A_need, rad_use=A_need / A_RAD[ct],
                T_sink=T_sink, T_cooler_in=c['T_cool_in'], T_comp_in=T1, running=True,
                part_load=load < 1 - 1e-9, load=load,
                m_gas=float(lp['m_g']), W_blower=float(lp['W_bl']), W_pump=float(lp['W_pump']),
                m_water=float(lp['m_w']), T_gas_hot=L['T_ex'], T_gas_cold=float(lp['T_g_out']),
                dp_loop_t=float(lp['dp']) / P_loop, K_gas=float(lp['m_g'])**2 / float(lp['rho_g']),
                Q_rec_max=float(lp['Q_rec_max']),
                T_rec=float(lp['T_rec']), T_rec_out=float(lp['T_ro']) if L['Q_rec'] > 0 else np.nan,
                T_w_hot=T_wh, T_w_cold=T_wc, **_hab_fields(hb, Q_hab))

def _plant_at(TIT, T_TES, t, ct, hb, T_sink, L):
    """Plant at a given turbine inlet temperature (coldest workable compressor inlet)."""
    r = cycle_row(TIT, ct)
    if not r['ok'].any():                                    # TIT too low for the cycle
        p = plant_off(T_TES, t, ct, L)
        p['too_cold'] = True
        return p
    Qh = hb['Q_heat']
    with np.errstate(invalid='ignore', divide='ignore'):
        lp_v = _loop(r, L)
        W_net_v = lp_v['W_net']
        A_need_v = _area_needed(r, lp_v['m_w'], T_sink, Qh)
        good = r['ok'] & (W_net_v > 0) & (A_need_v <= A_RAD[ct])
    if not good.any():
        # radiator cannot take full power at any T1: run at reduced power.
        # Waste heat (and area needed) scale with the cycle flow, so at each T1 the
        # radiator allows a fraction s = A_RAD / A_need of the full flow.
        with np.errstate(invalid='ignore', divide='ignore'):
            s = np.where(r['ok'] & np.isfinite(A_need_v) & (A_need_v > 0),
                         np.minimum(A_RAD[ct] / A_need_v, 1.0), 0.0)
            Wn = s * W_net_v
            lp_s = _loop_at_flow(r, Wn * 1e3 / r['w_net'], L)
            We = eta_gen * Wn - (W_aux_fixed + L['W_extra'] + lp_s['W_bl'] + lp_s['W_pump'])
        We = np.where(np.isfinite(We), We, -np.inf)
        if not (We > 0).any():                               # no net power at all: off
            return plant_off(T_TES, t, ct, L)
        # use the coldest T1 that still gives net power: best efficiency, so the
        # TES is drained as little as possible for the power the radiator allows
        i = int(np.argmax(We > 0))
        c = {k: float(r[k][i]) for k in KEYS}
        return _plant_state(T_TES, TIT, c, T1_grid[i], float(Wn[i]), A_RAD[ct],
                            T_sink, ct, load=float(s[i]), hb=hb, L=L)
    i = int(np.argmax(good))                                 # coldest T1 that works
    if i == 0:
        T1 = T1_grid[0]                                      # radiator has spare: bypass
    else:
        # interpolate where the radiator is exactly used up. Use 1/A_need, which
        # goes smoothly to 0 as T1 approaches the sink (A_need -> infinity there)
        with np.errstate(divide='ignore'):
            g0, g1 = 1 / A_need_v[i - 1], 1 / A_need_v[i]
        g = 1 / A_RAD[ct]
        if r['ok'][i - 1] and g0 < g <= g1:
            lo, hi = T1_grid[i - 1], T1_grid[i]                  # bisection: A_need(T1) = A_RAD
            for _ in range(20):
                mid = 0.5 * (lo + hi)
                cm = {k: float(np.interp(mid, T1_grid, r[k])) for k in KEYS}
                lm = _loop(cm, L)
                lo, hi = (lo, mid) if _area_needed(cm, lm['m_w'], T_sink, Qh) <= A_RAD[ct] else (mid, hi)
            T1 = hi
        else:
            T1 = T1_grid[i]
    c = {k: float(np.interp(T1, T1_grid, r[k])) for k in KEYS}
    lp = _loop(c, L)
    W_net = float(lp['W_net'])
    A_need = float(_area_needed(c, float(lp['m_w']), T_sink, Qh))
    return _plant_state(T_TES, TIT, c, T1, W_net, A_need, T_sink, ct, load=1.0, hb=hb, L=L)

def plant_point(T_TES, t=None, ct='simple', W_extra=0.0, Q_rec=0.0):
    """Plant at TES temperature T_TES. W_extra = heliostat drives (kWe, Step 9);
    Q_rec = solar heat into the gas loop this hour (kW, Step 9; 0 at night)."""
    T_sink = T_sink_night if t is None else float(T_sink_t(t))
    hb = habitat_balance(t)
    L = _loop_context(T_TES, W_extra, Q_rec, ct)
    TIT = min(L['T_ex'] - dT_PHX, TIT_top)
    return _plant_at(TIT, T_TES, t, ct, hb, T_sink, L)

def phx_UA(p, ct):
    """PHX size UA (W/K) needed at plant state p: UA = m_cycle * integral dh_c / (T_gas - T_cycle)."""
    if not p['running']:
        return np.nan
    h_c_out = np.interp(p['TIT'], _T_phx, _h_phx[ct])
    hc = p['h_phx_in'] + _f * (h_c_out - p['h_phx_in'])
    Tc = np.interp(hc, _h_phx[ct], _T_phx)
    hg_in = np.interp(p['T_gas_hot'], _T_g, _h_g)
    Tg = np.interp(hg_in - p['m_dot'] * (h_c_out - hc) / p['m_gas'], _h_g, _T_g)
    y = 1.0 / (Tg - Tc)
    return p['m_dot'] * float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(hc)))

# --- Lowest TES temperature at which the plant still makes full power (night) ---
for ct in cycles:
    T_ok = T_TES_max
    for T in np.arange(T_TES_max, T_TES_min - 0.1, -2.0):
        p = plant_point(T, None, ct)
        if not p['running'] or p['part_load']:
            break
        T_ok = T
    CYC[ct]['T_TES_cutoff_cycle'] = CYC[ct]['T_TES_cutoff']
    CYC[ct]['T_TES_cutoff'] = max(CYC[ct]['T_TES_cutoff'], T_ok)
    say(f"{ct:>11}: usable TES down to {CYC[ct]['T_TES_cutoff']:.0f} K "
        f"(turbine inlet = TES - {dT_TES_HX + dT_PHX:.0f} K at night)")

# --- Size the fixed radiator (night sink, coldest inlet, worst TES temperature) and the PHX ---
A_RAD_night_need, UA_PHX = {}, {}
for ct in cycles:
    pts = [plant_point(T, None, ct) for T in np.linspace(CYC[ct]['T_TES_cutoff'], T_TES_max, 30)[::-1]]
    need = [p['A_need'] for p in pts if np.isfinite(p['A_need'])]   # points where the plant runs
    A_RAD_night_need[ct] = max(need)
    A_RAD[ct] = rad_area_margin * A_RAD_night_need[ct]
    UA_PHX[ct] = max(phx_UA(p, ct) for p in pts if p['running'])
    p0 = plant_point(T_TES_max, None, ct)
    say(f"{ct:>11}: radiator for the TES floor {A_RAD[ct]:6.0f} m2 of panel, {rad_sides}-sided "
        f"({A_RAD[ct] * rad_areal_mass / 1e3:.1f} t at {rad_areal_mass:.0f} kg/m2) | PHX UA "
        f"{UA_PHX[ct]/1e3:.1f} kW/K | night, TES full: gas {p0['m_gas']:.2f} kg/s "
        f"{p0['T_gas_hot']:.0f}->{p0['T_gas_cold']:.0f} K; water {p0['m_water']:.2f} kg/s "
        f"{p0['T_w_cold']:.0f}->{p0['T_w_hot']:.0f} K (blower and radiator final in Step 10)")
