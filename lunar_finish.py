# lunar_finish.py -- finish the sizing for every cycle (Step 10c)
# Loaded by Step 10c of the notebook (keep this file next to the notebook), after Step 10b.
# Not a stand-alone script: it uses R (Step 10b), the sizing helpers (lunar_sizing.py,
# Step 10a), simulate() (Step 9), plant_point() and the radiator (Step 6) and the inputs.
# ---------------------------------------------------------------------
# For each cycle:
#   1. radiator resized to the peak of normal operation (x rad_area_margin)
#   2. final month run again with that radiator (from the sunrise state of Step 10b)
#   3. lowest end-A gas temperature that still gives full power at night
#   4. two sizing lines for the log
# ---------------------------------------------------------------------
def _finish(ct):
    n_log = len(RUN_LOG)
    r = R[ct]
    cy = r['cycle']
    _on = cy['running'] & ~cy['part_load']
    A_cons = rad_area_margin * A_RAD_night_need[ct]      # radiator from Step 6 (store floor)
    A_RAD[ct] = rad_area_margin * float(np.nanmax(cy['A_need'][_on])) * (1 + 1e-9)
    ok_r, cy = simulate(r['m_TES'], r['A_col'], r['T_sunrise'], n_cycle, t0=t_rise, ct=ct)
    cy.pop('T_end')
    r.update(cycle=cy, day={k: v[:n_day + 1] for k, v in cy.items()},
             night={k: v[n_day:] for k, v in cy.items()}, survived=r['survived'] and ok_r,
             periodic_error=cy['E_TES'][-1] - cy['E_TES'][0])
    # lowest gas temperature from end A at which the plant still gives full power at night
    T_ok = T_TES_max
    for T in np.arange(T_TES_max, T_TES_min - 0.1, -2.0):
        p = plant_point(T, None, ct)
        if not p['running'] or p['part_load']:
            break
        T_ok = T
    E_sol = np.sum(cy['Q_solar'][:-1]) * dt
    E_cur = np.sum(cy['Q_curtail'][:-1]) * dt
    say(f"{ct:>11}: TES {r['m_TES']/1e3:7.1f} t ({r['V_TES']:5.0f} m3, {N_segments} segments) | "
        f"collector {r['A_col']:5.0f} m2 | curtailed {100*E_cur/max(E_sol, 1e-9):4.1f}% | "
        f"survives {r['survived']} | repeats to {abs(r['periodic_error']):.1f} kWh")
    say(f"{'':>11}  blower 1 peak {np.nanmax(cy['W_bl1']):.1f} kWe, blower 2 peak "
        f"{np.nanmax(cy['W_blower']):.1f} kWe | radiator {A_cons:.0f} -> {A_RAD[ct]:.0f} m2, "
        f"peak use {np.nanmax(np.where(cy['running'], cy['rad_use'], np.nan))*100:.1f} % | "
        f"full power at night down to end-A gas {T_ok:.0f} K")
    log = RUN_LOG[n_log:]
    del RUN_LOG[n_log:]                                  # added back in cycle order below
    return dict(ct=ct, R=r, A_RAD=A_RAD[ct], A_RAD_conservative=A_cons, T_full=T_ok, log=log)

_radiator_from_step6()
_t0 = _time.time()
_res, _how = _run_part(_finish)
A_RAD_conservative = {}
for _r in _res:                                          # results back, in cycle order
    ct = _r['ct']
    R[ct] = _r['R']
    A_RAD[ct] = _r['A_RAD']
    A_RAD_conservative[ct] = _r['A_RAD_conservative']
    CYC[ct]['T_TES_full_power'] = _r['T_full']
    RUN_LOG.extend(_r['log'])
_part_time('10c (finish)', _t0, _how)
