# lunar_simulate.py -- hour-by-hour TES energy balance over the lunar day and night
# Loaded by Step 9 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses the regolith (Step 3), plant (Step 6) and solar
# (Step 8) functions, the inputs (Step 2) and np from Step 1.
# ---------------------------------------------------------------------
# TES energy balance over time (day and night), for cycle ct
#   dE_TES/dt = Q_solar(t, T_TES) - Q_curtail - Q_in(T_TES) - Q_loss(T_TES)
#   Q_solar   : solar heat into the TES gas loop, net of receiver losses (Step 8); the
#               receiver temperature comes from the gas loop (Step 6), so it is iterated
#   (the gas carries Q_solar from the receiver; it gives Q_in to the PHX and the rest to
#    the TES, so the TES change is Q_solar - Q_curtail - Q_in - Q_loss as before)
#   Q_curtail : solar heat that does not fit when the TES is full, or above what the gas loop
#               can carry at its design flow (heliostats defocused)
#   Q_in      : heat drawn by the cycle, including the extra heat for the parasitic loads
#               (Step 6); reduced in hours when the radiator cannot take full power
#   Q_loss    : heat lost through the insulation (Step 3)
# Parasitic loads: heliostat drives by day (w_helio_drive x mirror area); the gas-loop
# blower and the coolant pump day and night (Step 6). By day the blower moves the gas for
# the heat actually stored (Q_solar - Q_curtail - Q_in).
# Explicit time steps of dt hours. The arrays hold the state at the START of each step
# (t0, t0 + dt, ..., t0 + n_steps*dt: n_steps + 1 points) and the flows during that step;
# e_end is the TES energy after the last of these n_steps + 1 steps.
# Energy check: m_TES * (e_end - e_start) = sum over steps of
#   (Q_solar - Q_curtail - Q_in - Q_loss) * dt   (include Q_loss in any balance, Step 12)
# t0 = hours since sunrise at the first step
# require_full = True: the run also fails in any hour the plant cannot give the full
# 100 kWe (radiator limit, or blower too large near the TES floor); used for sizing
# ---------------------------------------------------------------------
e_full = e_regolith(T_TES_max)       # J/kg at full charge

def simulate(m_TES, A_col, e_start, n_steps, t0=0.0, ct='simple', require_full=True):
    keys = ['t', 'is_day', 'T_TES', 'TIT', 'eta', 'm_dot', 'Q_in', 'Q_out', 'Q_recup',
            'Q_hab', 'Q_rad', 'A_rad', 'T_sink', 'W_elec', 'E_TES', 'Q_solar', 'Q_curtail',
            'Q_sun_in', 'Q_rec_loss', 'eta_rec_t', 'W_net', 'W_gross', 'W_aux', 'Q_loss',
            'T_comp_in', 'A_need', 'rad_use', 'running', 'part_load', 'load', 'W_unserved',
            'Q_hab_need', 'Q_hab_unmet', 'Q_hab_cool', 'W_hp', 'hp_on', 'Q_hab_wall', 'Q_hab_rad', 'A_hab_rad',
            'dT_PHX_t', 'W_helio', 'W_blower', 'W_pump', 'm_gas', 'm_water', 'T_gas_hot',
            'T_gas_cold', 'T_rec', 'T_rec_out', 'T_w_hot', 'T_w_cold', 'dp_loop_t', 'K_gas']
    out = {k: [] for k in keys}
    T_cut = CYC[ct]['T_TES_cutoff']
    e = e_start
    survived = True
    T_rec = None                                             # receiver temperature, warm start
    for k in range(n_steps + 1):
        t = t0 + k * dt
        T = T_from_e(e)
        if T < T_cut - 1e-6:                                 # TES below the usable floor
            survived = False
            break
        Q_sun_in = A_col * float(solar_intercepted(t)) / 1e3   # kW concentrated light on receiver
        W_helio = w_helio_drive * A_col / 1e3 if Q_sun_in > 0 else 0.0   # kWe heliostat drives
        Q_loss = float(TES_heat_loss(T, m_TES))                # kW through the insulation
        room = (e_full - e) * m_TES / (1e3 * dt * 3600)        # kW the TES can still take this step
        # Solar heat kept, receiver temperature (from the gas loop), blower power and cycle
        # draw depend on each other: a few passes. When the TES is full the surplus is
        # curtailed (heliostats defocused) and the gas flow drops with it.
        if T_rec is None or not np.isfinite(T_rec):
            T_rec = T + dT_TES_HX + 0.5 * dT_charge + dT_wall    # first guess
        Q_in_g, Q_max_g = 0.0, np.inf
        for _ in range(12):
            Q_sun = A_col * float(solar_flux(t, T_rec)) / 1e3  # kW net into the gas loop
            Q_store = min(Q_sun, room + Q_in_g + Q_loss, Q_max_g)   # what fits (TES and gas loop)
            p = plant_point(T, t, ct, W_extra=W_helio, Q_rec=Q_store)
            if p['too_cold']:
                break
            T_new = p['T_rec'] if np.isfinite(p['T_rec']) else T_rec
            done = (abs(p['Q_in'] - Q_in_g) < 1e-4 and abs(T_new - T_rec) < 0.01
                    and abs(min(p['Q_rec_max'], 1e12) - min(Q_max_g, 1e12)) < 1e-3)
            Q_in_g, T_rec, Q_max_g = p['Q_in'], T_new, p['Q_rec_max']
            if done:
                break
        if p['too_cold']:                                    # turbine inlet below the cycle cutoff
            survived = False
            break
        if require_full and (not p['running'] or p['part_load']):   # demand not met this hour
            survived = False
            break
        Q_store = min(Q_sun, room + p['Q_in'] + Q_loss, p['Q_rec_max'])
        Q_curt = Q_sun - Q_store                               # kW curtailed (TES full or gas-loop flow limit)
        e_new = min(e + (Q_store - p['Q_in'] - Q_loss) * 1e3 * dt * 3600 / m_TES, e_full)
        out['t'].append(t)
        out['is_day'].append(np.mod(t, t_lunar) < t_day)
        out['E_TES'].append(m_TES * e / 3.6e6)               # kWh above T_TES_min
        out['Q_solar'].append(Q_sun)
        out['Q_curtail'].append(Q_curt)
        out['Q_loss'].append(Q_loss)
        out['W_helio'].append(W_helio)
        out['Q_sun_in'].append(Q_sun_in)
        out['Q_rec_loss'].append(Q_sun_in - Q_sun)            # reflection + radiation
        out['eta_rec_t'].append(Q_sun / Q_sun_in if Q_sun_in > 0 else np.nan)
        for key in ['T_TES', 'TIT', 'eta', 'm_dot', 'Q_in', 'Q_out', 'Q_recup',
                    'Q_hab', 'Q_rad', 'A_rad', 'T_sink', 'W_elec', 'W_net', 'W_gross', 'W_aux',
                    'T_comp_in', 'A_need', 'rad_use', 'running', 'part_load', 'load',
                    'Q_hab_need', 'Q_hab_unmet', 'Q_hab_cool', 'W_hp', 'hp_on', 'Q_hab_wall',
                    'Q_hab_rad', 'A_hab_rad', 'dT_PHX_t', 'W_blower', 'W_pump', 'm_gas',
                    'm_water', 'T_gas_hot', 'T_gas_cold', 'T_rec', 'T_rec_out', 'T_w_hot', 'T_w_cold',
                    'dp_loop_t', 'K_gas']:
            out[key].append(p[key])
        out['W_unserved'].append(max(W_elec_target - p['W_elec'], 0.0))   # kWe not delivered
        e = e_new
    out = {k: np.array(v) for k, v in out.items()}
    out['e_end'] = e
    return survived, out
