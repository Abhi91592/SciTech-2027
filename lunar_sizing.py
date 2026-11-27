# lunar_sizing.py -- sizing functions for the TES mass and the collector area (Step 10)
# Loaded by Step 10a of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses simulate() (Step 9), plant_point() and the radiator
# (Step 6), the sun timing (Step 8), the regolith energy (Step 3) and the inputs (Step 2).
# Defines only settings and functions; nothing runs until Steps 10b and 10c.
# ---------------------------------------------------------------------
# Sizing for each cycle (Eq. 21-25 + collector area, Section VI.3). Results in R[ct].
#   (a) Paper method (Eq. 21-24): heat at the design point x longest darkness, over the full
#       brick window T_TES_max -> T_TES_min (for comparison only)
#   (b) Smallest TES mass for the longest darkness alone: every segment at T_TES_max when it
#       starts (t_dark0, Step 8), no sun; the gas from end A must stay above the plant's
#       floor (CYC[ct]['T_TES_cutoff']) to its end
#   (c) Design TES mass = TES_mass_margin x (b), for dawn and dusk
#   (d) Smallest collector area that keeps the plant running month after month: start full
#       at t_dark0, run n_cycles_periodic design months with the real sun; the store must
#       never run empty, and its stored energy at the last darkness start must not be lower
#       than at the one before
#   (e) Final run: settle to the repeating state (ends at sunrise, t_rise), then one month
#       from sunrise. day = sunrise -> darkness start, night = the longest darkness.
# simulate(..., n_steps) runs n_steps + 1 steps (Step 9): a period of n steps uses n - 1.
# The radiator is the one Step 6 sized for the store floor; Step 10c resizes it afterwards.
# Bisections stop at a relative tolerance of size_rtol.
# ---------------------------------------------------------------------
size_rtol = 1e-2                       # relative tolerance of the mass and area bisections

def _bisect_min(works, lo, hi):
    """Smallest x in [lo, hi] with works(x) True (False below the answer, True above)."""
    while hi - lo > size_rtol * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if works(mid) else (mid, hi)
    return hi

def size_system(ct):
    """Size the plant for cycle ct: TES mass, collector area and the final month."""
    # (a) Paper method: design-point Q_in for the whole darkness, full brick window
    p_design = plant_point(T_TES_max, ct=ct)
    E_req_design = p_design['Q_in'] * t_night                         # kWh, Eq. 21
    dE_window = e_regolith(T_TES_max) / 3.6e6                         # kWh/kg (relative to T_TES_min)
    m_TES_design = E_req_design / dE_window                           # Eq. 24

    # (b) Smallest TES mass for the longest darkness alone (start full, no sun)
    def survives_night(m):
        return simulate(m, 0.0, T_TES_max, n_dark - 1, t0=t_dark0, ct=ct)[0]
    lo, hi = 0.1 * m_TES_design, 20 * m_TES_design
    if survives_night(lo) or not survives_night(hi):
        raise RuntimeError(f"{ct}: night-mass search range does not bracket the answer")
    m_TES_night = _bisect_min(survives_night, lo, hi)

    # (c) Design TES mass with margin
    m_TES = TES_mass_margin * m_TES_night
    V_TES = m_TES / rho_TES                                           # Eq. 25

    # (d) Smallest collector area that does not run the store down month to month
    N = n_cycles_periodic
    def survives_periodic(A):
        ok, out = simulate(m_TES, A, T_TES_max, N * n_cycle, t0=t_dark0, ct=ct)
        if not ok:
            return False
        E = out['E_TES']                                              # index j*n_cycle = j-th darkness start
        return E[N * n_cycle] >= E[(N - 1) * n_cycle] - 1e-6
    A_lo, A_hi = 0.0, 1e3
    while not survives_periodic(A_hi):
        A_lo, A_hi = A_hi, 2 * A_hi
        if A_hi > 1e8:
            raise RuntimeError(f"{ct}: no collector area works, increase TES_mass_margin")
    A_col = _bisect_min(survives_periodic, A_lo, A_hi)

    # (e) Final run: settle (ends at sunrise, t_rise), then one month from sunrise
    ok_settle, pre = simulate(m_TES, A_col, T_TES_max, (N - 1) * n_cycle + n_dark - 1,
                              t0=t_dark0, ct=ct)
    T_sunrise = pre['T_end']                                          # segment temperatures at sunrise
    survived, cycle = simulate(m_TES, A_col, T_sunrise, n_cycle, t0=t_rise, ct=ct)
    cycle.pop('T_end')
    day = {k: v[:n_day + 1] for k, v in cycle.items()}                # sunrise -> darkness start
    night = {k: v[n_day:] for k, v in cycle.items()}                  # longest darkness
    return dict(p_design=p_design, E_req_design=E_req_design, m_TES_design=m_TES_design,
                T_sunrise=T_sunrise, m_TES_night=m_TES_night, m_TES=m_TES, V_TES=V_TES,
                A_col=A_col, cycle=cycle, day=day, night=night,
                survived=survived and ok_settle,
                periodic_error=cycle['E_TES'][-1] - cycle['E_TES'][0])   # kWh, one month

# --- Run one sizing part for every cycle ----------------------------------------------
# The cycles are independent, so with PARALLEL = True they are sized at the same time in
# separate processes ('fork': Linux, Colab, macOS). On Windows, or with PARALLEL = False,
# the cycles are sized one after the other. Same results either way.
import multiprocessing as _mp, warnings as _warnings, time as _time
PARALLEL = True

def _run_part(fn):
    """fn(ct) for every cycle, in cycle order: (results, how it ran)."""
    if PARALLEL and len(cycles) > 1 and 'fork' in _mp.get_all_start_methods():
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore', DeprecationWarning)   # fork from a threaded kernel
            with _mp.get_context('fork').Pool(len(cycles)) as pool:
                return pool.map(fn, cycles), 'both cycles in parallel'
    return [fn(ct) for ct in cycles], 'one cycle after the other'

def _part_time(label, t0, how):
    msg = f"Step {label}: {(_time.time() - t0)/60:.1f} min ({how})"
    say(msg)
    print(msg)

def _radiator_from_step6():
    """Radiator as sized in Step 6 (for the store floor): the sizing uses this one."""
    for ct in cycles:
        A_RAD[ct] = rad_area_margin * A_RAD_night_need[ct]
