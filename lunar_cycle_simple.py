# lunar_cycle_simple.py -- simple sCO2 Brayton cycle and its (T1, TIT) table
# Loaded by Step 4 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses np, CP, say() (Step 1) and the inputs (Step 2).
# Defines: P_simple, P_phx_out, compressor(), turbine(), simple_cycle(), T1_grid, TIT_grid,
#          KEYS, build_cycle_table(), cycle_row(), cycle_at(), CYC (Step 5 adds 'recuperated')
# ---------------------------------------------------------------------
# Simple sCO2 Brayton cycle (Appendix Table A, Eq. 1-10) with pressure drops
# States: 1 compressor in, 2 compressor out / PHX in,
#         3 PHX out / turbine in, 4 turbine out / cooler in
# Pressures: P1 = P_low, P2 = P_high, P3 = P2*(1 - dp_PHX),
#            P4 = P1/(1 - dp_cooler)  (the turbine expands only to the pressure
#            that leaves P_low after the cooler/radiator)
# Isentropic efficiencies are applied to ENTHALPY (not temperature, which is only
# valid for an ideal gas with constant cp); outlet temperatures follow from (h, P):
#   compressor: h2 = h1 + (h2s - h1) / eta_comp,   T2 = T(h2, P2)      Eq. 8
#   turbine:    h4 = h3 - eta_turb * (h3 - h4s),   T4 = T(h4, P4)      Eq. 10
# T1 = compressor inlet temperature, set each hour by the radiator (Step 6).
# Results are pre-computed on a (T1, TIT) grid and stored in CYC[cycle];
# cycle_row(TIT, ct) gives all T1 values at one TIT, cycle_at(TIT, ct, T1) one point.
# Validated against Sandia SAND2012-6426 (see the validation note in the project).
# ---------------------------------------------------------------------
P_simple = dict(P1=P_low, P2=P_high, P3=P_high * (1 - dp_PHX), P4=P_low / (1 - dp_cooler))
P_phx_out = {'simple': P_simple['P3']}               # PHX outlet = turbine inlet pressure

# Speed: the compressor depends only on T1 and the turbine only on TIT, so each is
# calculated once per value and reused across the grid (same CoolProp calls, same
# numbers: the results are identical). The cache keys include every input used, so a
# changed pressure, efficiency or fluid is recalculated.
from functools import lru_cache

@lru_cache(maxsize=None)
def compressor(T1, P_in, P_out, eta, fluid):
    """Compressor (Eq. 8): isentropic efficiency on enthalpy -> (h1, h2, T2)."""
    h1 = CP.PropsSI('H', 'T', T1, 'P', P_in, fluid)
    s1 = CP.PropsSI('S', 'T', T1, 'P', P_in, fluid)
    h2s = CP.PropsSI('H', 'S', s1, 'P', P_out, fluid)
    h2 = h1 + (h2s - h1) / eta
    return h1, h2, CP.PropsSI('T', 'H', h2, 'P', P_out, fluid)

@lru_cache(maxsize=None)
def turbine(TIT, P_in, P_out, eta, fluid):
    """Turbine (Eq. 10): isentropic efficiency on enthalpy -> (h_in, h_out, T_out)."""
    h3 = CP.PropsSI('H', 'T', TIT, 'P', P_in, fluid)
    s3 = CP.PropsSI('S', 'T', TIT, 'P', P_in, fluid)
    h4s = CP.PropsSI('H', 'S', s3, 'P', P_out, fluid)
    h4 = h3 - eta * (h3 - h4s)
    return h3, h4, CP.PropsSI('T', 'H', h4, 'P', P_out, fluid)

def simple_cycle(TIT, T1=None):
    T1 = T_comp_in if T1 is None else T1
    P = P_simple
    # Compressor (Eq. 8) and turbine (Eq. 10): isentropic efficiencies on enthalpy
    h1, h2, T2 = compressor(float(T1), P['P1'], P['P2'], eta_comp, FLUID)
    h3, h4, T4 = turbine(float(TIT), P['P3'], P['P4'], eta_turb, FLUID)

    w_c = h2 - h1                 # Eq. 1 (per kg)
    q_in = h3 - h2                # Eq. 2
    w_t = h3 - h4                 # Eq. 3
    q_out = h4 - h1               # Eq. 4
    w_net = w_t - w_c             # Eq. 5
    eta = w_net / q_in            # Eq. 6
    return dict(eta=eta, w_net=w_net, w_c=w_c, w_t=w_t, q_in=q_in, q_out=q_out,
                q_rec=0.0, T1=T1, T2=T2, T4=T4, T_cool_in=T4, h_cool_in=h4, h_phx_in=h2,
                h1=h1, h2=h2, h3=h3, h4=h4,
                P=[P['P1'], P['P2'], P['P3'], P['P4']])

# --- Grids: compressor inlet T1 x turbine inlet TIT ---
T1_grid = np.arange(T_comp_in, T_comp_in_max + 0.1, 2.5)
TIT_grid = np.unique(np.concatenate([np.arange(T_comp_in + 30, TIT_top, 5.0), [TIT_max, TIT_top]]))
KEYS = ['eta', 'w_net', 'w_c', 'w_t', 'q_in', 'q_out', 'q_rec', 'h1', 'h_cool_in', 'T_cool_in',
        'h_phx_in']

def build_cycle_table(cycle_fn):
    """Pre-compute a cycle on the (T1, TIT) grid. The TES cutoff uses the coldest T1."""
    G = {k: np.full((len(T1_grid), len(TIT_grid)), np.nan) for k in KEYS}
    for i, T1 in enumerate(T1_grid):
        for j, TIT in enumerate(TIT_grid):
            try:
                r = cycle_fn(TIT, T1)
            except ValueError:            # CoolProp could not solve this point: left as NaN
                continue
            for k in KEYS:
                G[k][i, j] = r[k]
    with np.errstate(invalid='ignore'):
        ok = (G['w_net'][0] > 0) & (G['q_in'][0] > 0) & (G['eta'][0] >= eta_cycle_min)
    if not ok[-1]:
        raise RuntimeError(f"{cycle_fn.__name__}: the cycle does not work even at the top of the "
                           f"TIT grid ({TIT_grid[-1]:.0f} K); check T_TES_max, pressures and efficiencies")
    # lowest TIT at which the cycle still works (scan down from the top, stop at first failure)
    first_fail_from_top = np.argmax(~ok[::-1])
    TIT_cut = TIT_grid[len(TIT_grid) - first_fail_from_top] if not ok.all() else TIT_grid[0]
    return dict(G=G, TIT_cutoff=TIT_cut, T_TES_cutoff=max(T_TES_min, TIT_cut + dT_PHX))   # gas leaves the TES at T_TES

def cycle_row(TIT, ct):
    """All cycle results versus T1 (arrays over T1_grid) at one TIT, plus a 'works' mask.
    TIT outside the grid is held at the nearest grid edge (no extrapolation)."""
    TIT = float(np.clip(TIT, TIT_grid[0], TIT_grid[-1]))
    j = int(np.clip(np.searchsorted(TIT_grid, TIT) - 1, 0, len(TIT_grid) - 2))
    w = (TIT - TIT_grid[j]) / (TIT_grid[j + 1] - TIT_grid[j])
    G = CYC[ct]['G']
    r = {k: G[k][:, j] * (1 - w) + G[k][:, j + 1] * w for k in KEYS}
    with np.errstate(invalid='ignore'):
        r['ok'] = (r['w_net'] > 0) & (r['q_in'] > 0) & (r['eta'] >= eta_cycle_min)
    return r

def cycle_at(TIT, ct='simple', T1=None):
    """Cycle results at one (TIT, T1); T1 defaults to the coldest compressor inlet."""
    T1 = T_comp_in if T1 is None else T1
    r = cycle_row(TIT, ct)
    return {k: float(np.interp(T1, T1_grid, r[k])) for k in KEYS}

say(f"Building simple-cycle map ({len(T1_grid)} x {len(TIT_grid)} points)...")
CYC = {'simple': build_cycle_table(simple_cycle)}
say(f"simple cycle: turbine {P_simple['P3']/1e6:.2f} -> {P_simple['P4']/1e6:.2f} MPa, "
    f"eta at TIT {TIT_max:.0f} K / T1 {T_comp_in:.0f} K = {cycle_at(TIT_max, 'simple')['eta']*100:.2f} %, "
    f"usable TES down to {CYC['simple']['T_TES_cutoff']:.0f} K")