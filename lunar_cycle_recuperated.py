# lunar_cycle_recuperated.py -- recuperated sCO2 Brayton cycle and its (T1, TIT) table
# Loaded by Step 5 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses np, CP, say() (Step 1), the inputs (Step 2) and the
# compressor/turbine functions, grids and CYC from Step 4 (lunar_cycle_simple.py).
# Defines: P_recup, recup_min_dT(), recuperated_cycle(), cycle_state_table(),
#          CYC['recuperated']
# ---------------------------------------------------------------------
# Recuperated sCO2 Brayton cycle (Appendix Table A, Eq. 11-20)
# States: 1 compressor in        2 compressor out / recuperator cold in
#         3 recuperator cold out / PHX in      4 PHX out / turbine in (TIT)
#         5 turbine out / recuperator hot in   6 recuperator hot out / cooler in
# Pressures: P1 = P_low, P2 = P_high, P3 = P2(1-dp_recup), P4 = P3(1-dp_PHX),
#            P6 = P1/(1-dp_cooler), P5 = P6/(1-dp_recup)
# Compressor and turbine as in Step 4 (isentropic efficiencies on ENTHALPY, Eq. 8, 10).
# Recuperator (Eq. 11-14), counterflow:
#   q_max = the largest duty for which the hot stream stays hotter than the cold stream
#           EVERYWHERE inside the exchanger, not only at its two ends. sCO2 cp changes
#           strongly with temperature, so the closest approach can sit inside the
#           exchanger; the profile is checked at n_recup_check points and q_max is
#           reduced (bisection) if the streams would cross.
#   q_rec = eps_recup x q_max;   h3 = h2 + q_rec,  h6 = h5 - q_rec   (Q_cold = Q_hot)
#   The recuperator is always in the flow path (no bypass). If the turbine exhaust (5) is
#   colder than the compressor outlet (2), heat flows the other way, from the high-pressure
#   stream to the exhaust: q_rec < 0 with the same effectiveness, so the sCO2 reaches the
#   PHX colder and the precooler hotter. The duty passes smoothly through 0 at T5 = T2.
# The smallest temperature difference inside the exchanger (between whichever stream is
# hotter and the other) is returned as dT_pinch (K) when pinch=True (report only).
# Validated against Sandia SAND2012-6426 (see the validation note in the project).
# ---------------------------------------------------------------------
P_recup = dict(P1=P_low, P2=P_high, P3=P_high * (1 - dp_recup))
P_recup['P4'] = P_recup['P3'] * (1 - dp_PHX)                 # turbine inlet
P_recup['P6'] = P_low / (1 - dp_cooler)                      # recuperator hot out / cooler in
P_recup['P5'] = P_recup['P6'] / (1 - dp_recup)               # turbine outlet
P_phx_out['recuperated'] = P_recup['P4']
n_recup_check = 10                                           # points along the exchanger

# Speed: the recuperator crossing check (called ~22,000 times, ~94 % of this cell's time)
# reads T(h, P) from CoolProp's property tables (BICUBIC&HEOS: built once, saved on disk,
# ~300x faster, T within 0.001 K of the full equation of state). Every cycle state
# (compressor, turbine, recuperator ends) still uses the full equation of state (PropsSI).
_AS_tab = CP.AbstractState('BICUBIC&HEOS', FLUID)

def _T_hP_tab(h, p):
    out = np.empty(len(h))
    for i in range(len(h)):
        _AS_tab.update(CP.HmassP_INPUTS, h[i], p[i])
        out[i] = _AS_tab.T()
    return out

def recup_min_dT(h2, h5, q, P):
    """Smallest temperature difference, in the direction the heat flows, along a counterflow
    recuperator carrying duty q (J/kg; q > 0: exhaust -> high-pressure side, q < 0: the
    reverse). Point f runs from f = 0 (exhaust outlet h5 - q, high-pressure inlet h2) to
    f = 1 (exhaust inlet h5, high-pressure outlet h2 + q); pressures vary linearly."""
    f = np.linspace(0.0, 1.0, n_recup_check + 1)
    T_lp = _T_hP_tab(h5 - q + f * q, P['P6'] + f * (P['P5'] - P['P6']))
    T_hp = _T_hP_tab(h2 + f * q, P['P2'] + f * (P['P3'] - P['P2']))
    sgn = 1.0 if q >= 0 else -1.0
    return float(np.min(sgn * (np.asarray(T_lp) - np.asarray(T_hp))))

def recuperated_cycle(TIT, T1=None, pinch=False):
    """pinch=True also returns dT_pinch (slower; used for the report, not for the tables)."""
    T1 = T_comp_in if T1 is None else T1
    P = P_recup
    # Compressor (Eq. 8) and turbine (Eq. 10): same functions as Step 4, each calculated
    # once per T1 / TIT and reused (identical numbers)
    h1, h2, T2 = compressor(float(T1), P['P1'], P['P2'], eta_comp, FLUID)
    h4, h5, T5 = turbine(float(TIT), P['P4'], P['P5'], eta_turb, FLUID)

    # Recuperator (Eq. 11-14)
    # Largest duty from the end temperatures: the stream giving heat can at most be cooled to
    # the other stream's inlet, the stream taking heat at most heated to the other's inlet.
    dh_lp = h5 - CP.PropsSI('H', 'T', T2, 'P', P['P6'], FLUID)        # exhaust brought to T2
    dh_hp = CP.PropsSI('H', 'T', T5, 'P', P['P3'], FLUID) - h2        # HP stream brought to T5
    if T5 >= T2:                                   # normal: exhaust heats the HP stream
        q_max = max(min(dh_lp, dh_hp), 0.0)
    else:                                          # reverse: HP stream heats the exhaust
        q_max = -max(min(-dh_lp, -dh_hp), 0.0)     # q_max <= 0
    if q_max != 0 and recup_min_dT(h2, h5, q_max, P) < 0:
        # the streams would cross inside the exchanger: largest |duty| with no crossing
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if recup_min_dT(h2, h5, mid * q_max, P) >= 0 else (lo, mid)
        q_max = lo * q_max
    q_rec = eps_recup * q_max                            # Eq. 13: Q_cold = Q_hot
    h3 = h2 + q_rec                                      # Eq. 11
    h6 = h5 - q_rec                                      # Eq. 12
    T3 = CP.PropsSI('T', 'H', h3, 'P', P['P3'], FLUID)
    T6 = CP.PropsSI('T', 'H', h6, 'P', P['P6'], FLUID)
    dT_pinch = recup_min_dT(h2, h5, q_rec, P) if (pinch and q_rec != 0) else np.nan

    w_c = h2 - h1                 # Eq. 17 (per kg)
    q_in = h4 - h3                # Eq. 15: heat from TES after recuperation
    w_t = h4 - h5                 # Eq. 16
    q_out = h6 - h1               # Eq. 20: heat rejected after recuperation
    w_net = w_t - w_c             # Eq. 18
    eta = w_net / q_in            # Eq. 19
    return dict(eta=eta, w_net=w_net, w_c=w_c, w_t=w_t, q_in=q_in, q_out=q_out,
                q_rec=q_rec, dT_pinch=dT_pinch,
                T1=T1, T2=T2, T3=T3, T5=T5, T6=T6, T_cool_in=T6, h_cool_in=h6, h_phx_in=h3,
                h1=h1, h2=h2, h3=h3, h4=h4, h5=h5, h6=h6,
                P=[P['P1'], P['P2'], P['P3'], P['P4'], P['P5'], P['P6']])

def cycle_state_table(TIT, ct, T1=None):
    """Cycle results and (name, T, P) of every state, for the report."""
    if ct == 'recuperated':
        r = recuperated_cycle(TIT, T1, pinch=True)
        names = ['1  Compressor inlet', '2  Compressor outlet', '3  Recuperator cold out',
                 '4  Turbine inlet', '5  Turbine outlet', '6  Cooler inlet']
        temps = [r['T1'], r['T2'], r['T3'], TIT, r['T5'], r['T6']]
    else:
        r = simple_cycle(TIT, T1)
        names = ['1  Compressor inlet', '2  Compressor outlet', '3  Turbine inlet', '4  Turbine outlet']
        temps = [r['T1'], r['T2'], TIT, r['T4']]
    return r, list(zip(names, temps, r['P']))

say(f"Building recuperated-cycle map ({len(T1_grid)} x {len(TIT_grid)} points)...")
CYC['recuperated'] = build_cycle_table(recuperated_cycle)
_r = recuperated_cycle(TIT_max, T_comp_in, pinch=True)
say(f"recuperated cycle: turbine {P_recup['P4']/1e6:.2f} -> {P_recup['P5']/1e6:.2f} MPa, "
    f"eta at TIT {TIT_max:.0f} K / T1 {T_comp_in:.0f} K = {cycle_at(TIT_max, 'recuperated')['eta']*100:.2f} %, "
    f"recuperator pinch {_r['dT_pinch']:.1f} K, "
    f"usable TES down to {CYC['recuperated']['T_TES_cutoff']:.0f} K")
