# lunar_discharge_table.py -- TES discharge table (night), both cycles
# Loaded by Step 7 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses plant_point() and the cycle maps from Steps 4-6.
# ---------------------------------------------------------------------
# Discrete TES discharge temperatures (Section VI.4), both cycles (night).
# TIT = T_TES - dT_TES_HX - dT_PHX at night (TES gas loop and PHX, Step 6).
#   eta_cyc = W_net / Q_in            (cycle: shaft power / heat from the TES)
#   eta_net = W_elec / Q_in           (plant: electricity to the habitat / heat from the TES,
#                                      after the generator and the parasitic loads)
#   gap     = T_TES - TIT             (TES heat exchanger + PHX approach)
# ---------------------------------------------------------------------
say(f"{W_elec_target:.0f} kWe to the habitat, PR = {P_high/P_low:.2f}, regolith = {regolith_type} (night)\n")
_w = 31
say(f"{'':6} {'---------- Simple ----------':>{_w}}   {'-------- Recuperated --------':>{_w}}")
say(f"{'T_TES':>6} {'eta_cyc':>8} {'eta_net':>8} {'Q_in':>7} {'gap':>6}   "
    f"{'eta_cyc':>8} {'eta_net':>8} {'Q_in':>7} {'gap':>6}")
say(f"{'K':>6} {'%':>8} {'%':>8} {'kW':>7} {'K':>6}   {'%':>8} {'%':>8} {'kW':>7} {'K':>6}")
_worse = []
for T in np.arange(T_TES_max, T_TES_min - 1, -25):
    cols, Q = [], {}
    for ct in cycles:
        p = plant_point(T, ct=ct) if T >= CYC[ct]['T_TES_cutoff'] else None
        if p is None or not p['running']:
            cols.append(f"{'cannot run':>{_w}}")
        else:
            Q[ct] = p['Q_in']
            cols.append(f"{p['eta']*100:8.1f} {p['W_elec']/p['Q_in']*100:8.1f} "
                        f"{p['Q_in']:7.0f} {p['dT_PHX_t']:6.1f}")
    if len(Q) == 2 and Q['recuperated'] > Q['simple']:
        _worse.append(T)
    say(f"{T:6.0f} " + "   ".join(cols))
for ct in cycles:
    say(f"Usable TES window ({ct}): {T_TES_max:.0f} K down to {CYC[ct]['T_TES_cutoff']:.0f} K")
if _worse:
    say(f"Note: at TES <= {max(_worse):.0f} K the turbine exhaust is close to or colder than the "
        f"compressor outlet, so the recuperator transfers little heat, or passes heat backwards to the "
        f"exhaust, and adds pressure drop; the recuperated cycle then draws MORE heat than the simple cycle.")
