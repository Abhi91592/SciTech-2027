# lunar_discharge_table.py -- discharge table (no sun), both plants
# Loaded by Step 7 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses plant_point() and the cycle maps from Steps 4-6.
# ---------------------------------------------------------------------
# Plant performance with no sun, every 25 K of the gas leaving TES end A (Section VI.4).
# TIT = T_gas_A - dT_PHX (discharging gas loop and PHX, Step 6).
#   eta_cyc = W_net / Q_in            (cycle: shaft power / heat from the TES)
#   eta_net = W_elec / Q_in           (plant: electricity to the habitat / heat from the TES,
#                                      after the generator and all parasitic loads)
#   blower2 = discharging blower power, kWe
#   T_B     = gas going back into TES end B (sets how cold the bricks there can get), K
# ---------------------------------------------------------------------
say(f"{W_elec_target:.0f} kWe to the habitat, PR = {P_high/P_low:.2f}, regolith = {regolith_type}, "
    f"TES loops {P_loop/1e6:.0f} MPa (no sun)\n")
_w = 38
say(f"{'':6} {'-------------- Simple --------------':>{_w}}   {'----------- Recuperated -----------':>{_w}}")
say(f"{'T_gas_A':>7} {'eta_cyc':>8} {'eta_net':>8} {'Q_in':>7} {'blower2':>8} {'T_B':>5}   "
    f"{'eta_cyc':>8} {'eta_net':>8} {'Q_in':>7} {'blower2':>8} {'T_B':>5}")
say(f"{'K':>7} {'%':>8} {'%':>8} {'kW':>7} {'kWe':>8} {'K':>5}   {'%':>8} {'%':>8} {'kW':>7} {'kWe':>8} {'K':>5}")
_worse = []
for T in np.arange(T_TES_max, T_TES_min - 1, -25):
    cols, Q = [], {}
    for ct in cycles:
        p = plant_point(T, ct=ct) if T >= CYC[ct]['T_TES_cutoff'] else None
        if p is None or not p['running']:
            cols.append(f"{'cannot run':>{_w}}")
        else:
            Q[ct] = p['Q_in']
            cols.append(f"{p['eta']*100:8.1f} {p['W_elec']/p['Q_in']*100:8.1f} {p['Q_in']:7.0f} "
                        f"{p['W_blower']:8.1f} {p['T_gas_cold']:5.0f}")
    if len(Q) == 2 and Q['recuperated'] > Q['simple']:
        _worse.append(T)
    say(f"{T:7.0f} " + "   ".join(cols))
for ct in cycles:
    say(f"Usable gas from TES end A ({ct}): {T_TES_max:.0f} K down to {CYC[ct]['T_TES_cutoff']:.0f} K")
if _worse:
    say(f"Note: with the gas from end A at or below {max(_worse):.0f} K the turbine exhaust is close to or "
        f"colder than the compressor outlet, so the recuperator transfers little heat, or passes heat "
        f"backwards to the exhaust, and adds pressure drop; the recuperated cycle then draws MORE heat "
        f"than the simple cycle.")