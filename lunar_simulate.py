# lunar_simulate.py -- step-by-step run of the plant over the lunar month (two TES gas loops)
# Loaded by Step 9 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses the TES segments (Step 3), the plant (Step 6), the sun
# (Step 8), the inputs (Step 2) and np from Step 1.
# ---------------------------------------------------------------------
# State: the brick temperatures of the N_segments TES segments (segment 0 = END A, hot;
# segment N-1 = END B, cold). The gas leaves each segment at that segment's temperature, so
#   gas from END A (to the PHX)      = T_seg[0]
#   gas from END B (to the receiver) = T_seg[-1]
# Each step t (the sun is taken at the middle of the step, t + dt/2):
#   1. TES empty?  gas from end A below the plant's floor (CYC[ct]['T_TES_cutoff']) -> fail
#   2. Receiver:   inlet = gas from end B, outlet = T_rec_set (held by blower 1),
#                  T_rec = (inlet + T_rec_set)/2 + dT_wall;  Q_solar = A_col * solar_flux() (kW)
#   3. Charging:   the store takes Q_solar up to what still fits this step
#                  ("full" = every segment at T_TES_max); the rest is curtailed (defocused).
#                  Blower 1 flow = Q_charge / (h(T_rec_set) - h(gas from end B)), capped at the
#                  flow that gives a dT_charge rise in the receiver; with a hotter end B the
#                  store takes less heat and the rest is curtailed.
#                  Extra loads for the plant: heliostat drives (sun visible) + blower 1.
#   4. Plant:      plant_point(gas from end A, t, ct, extra loads) -> heat to the cycle,
#                  blower 2 flow, gas back into end B. (Blower 1 and the plant depend on each
#                  other through the parasitic loads: a few passes.)
#   5. TES:        advance the segments over dt with both flows (lunar_tes.advance).
# Arrays hold the state at the START of each step (t0, t0 + dt, ..., n_steps + 1 points) and
# the flows during that step; T_end = segment temperatures after the last step.
# Energy check: stored energy change = sum over steps of (Q_ch - Q_dis - Q_loss) * dt
#   Q_ch = heat into the bricks from the charging gas, Q_dis = heat taken by the discharging
#   gas (= heat to the cycle, Q_in, up to the change of the end-A gas within a step).
# require_full = True: the run also fails in any step the plant cannot give the full 100 kWe.
# ---------------------------------------------------------------------
_PLANT_KEYS = ['TIT', 'eta', 'm_dot', 'Q_in', 'Q_out', 'Q_recup', 'Q_hab', 'Q_rad', 'A_rad',
               'T_sink', 'W_elec', 'W_net', 'W_gross', 'W_aux', 'T_comp_in', 'A_need', 'rad_use',
               'running', 'part_load', 'load', 'Q_hab_need', 'Q_hab_unmet', 'Q_hab_cool', 'W_hp',
               'hp_on', 'Q_hab_wall', 'Q_hab_rad', 'A_hab_rad', 'dT_PHX_t', 'W_blower', 'W_pump',
               'm_gas', 'm_water', 'T_gas_hot', 'T_gas_cold', 'T_w_hot', 'T_w_cold']

def simulate(m_TES, A_col, T_start, n_steps, t0=0.0, ct='simple', require_full=True):
    geo = tes_geometry(m_TES)
    T_s = np.broadcast_to(np.asarray(T_start, dtype=float), (N_segments,)).copy()
    E_full = stored_energy(np.full(N_segments, T_TES_max), geo)          # kWh
    T_cut = CYC[ct]['T_TES_cutoff']
    h_set = gas_h(T_rec_set)
    keys = ['t', 'is_day', 'T_seg', 'T_A', 'T_B', 'E_TES', 'Q_sun_in', 'Q_solar', 'Q_curtail',
            'Q_ch', 'Q_dis', 'Q_loss', 'Q_rec_loss', 'eta_rec_t', 'T_rec', 'm_ch', 'W_bl1',
            'W_helio', 'W_unserved'] + _PLANT_KEYS
    out = {k: [] for k in keys}
    survived = True
    for k in range(n_steps + 1):
        t = t0 + k * dt
        ts = t + 0.5 * dt                                     # sun at mid-step
        T_A, T_B = T_s[0], T_s[-1]
        if T_A < T_cut - 1e-6:                                # TES empty for this plant
            survived = False
            break
        # --- sun and receiver
        Q_sun_in = A_col * float(solar_intercepted(ts)) / 1e3          # kW light on receiver
        T_rec = 0.5 * (T_B + T_rec_set) + dT_wall
        Q_solar = A_col * float(solar_flux(ts, T_rec)) / 1e3           # kW net into the gas
        W_helio = w_helio_drive * A_col / 1e3 if Q_sun_in > 0 else 0.0
        E_now = stored_energy(T_s, geo)
        Q_loss_est = float(heat_loss(T_s, geo).sum()) / 1e3
        # --- charging flow and plant (parasitic loads couple them: a few passes)
        W_bl1, m_ch, Q_ch_set = 0.0, 0.0, 0.0
        for _ in range(6):
            p = plant_point(T_A, ts, ct, W_extra=W_helio + W_bl1)
            Q_room = (E_full - E_now) / dt + (p['Q_in'] if p['running'] else 0.0) + Q_loss_est
            Q_ch_set = max(0.0, min(Q_solar, Q_room))
            m_ch = (Q_ch_set * 1e3 / (h_set - gas_h(min(T_B, T_rec_set - dT_charge)))
                    if Q_ch_set > 0 and T_B < T_rec_set - 1e-6 else 0.0)
            W_new = blower_power(m_ch, T_B)
            if abs(W_new - W_bl1) < 1e-4:
                break
            W_bl1 = W_new
        if p['too_cold'] or (require_full and (not p['running'] or p['part_load'])):
            survived = False
            break
        # --- move the TES segments over the step
        m_dis = p['m_gas'] if p['running'] else 0.0
        T_new, o = advance(T_s, dt, geo, m_ch=m_ch, T_ch_in=T_rec_set,
                           m_dis=m_dis, T_dis_in=p['T_gas_cold'])
        # --- record
        out['t'].append(t)
        out['is_day'].append(bool(sun_visible(ts) > 0))
        out['T_seg'].append(T_s.copy())
        out['T_A'].append(T_A)
        out['T_B'].append(T_B)
        out['E_TES'].append(E_now)
        out['Q_sun_in'].append(Q_sun_in)
        out['Q_solar'].append(Q_solar)
        out['Q_ch'].append(o['Q_ch'])
        out['Q_curtail'].append(max(Q_solar - o['Q_ch'], 0.0))
        out['Q_dis'].append(o['Q_dis'])
        out['Q_loss'].append(o['Q_loss'])
        out['Q_rec_loss'].append(Q_sun_in - Q_solar)                   # reflection + radiation
        out['eta_rec_t'].append(Q_solar / Q_sun_in if Q_sun_in > 0 else np.nan)
        out['T_rec'].append(T_rec if Q_sun_in > 0 else np.nan)
        out['m_ch'].append(m_ch)
        out['W_bl1'].append(W_bl1)
        out['W_helio'].append(W_helio)
        out['W_unserved'].append(max(W_elec_target - p['W_elec'], 0.0))
        for key in _PLANT_KEYS:
            out[key].append(p[key])
        T_s = T_new
    out = {k: np.array(v) for k, v in out.items()}
    out['T_end'] = T_s
    return survived, out