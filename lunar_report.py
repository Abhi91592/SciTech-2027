# lunar_report.py -- report, figures, references and the results PDF, in one go
# Loaded by Step 12 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses the results of Steps 1-11 (R, CYC, A_RAD, ...).
#   Part A  text report                       (was notebook Step 12)
#   Part B  AIAA figures, saved to FIG_DIR    (was Step 13)
#   Part C  references, assumptions, validation (was Step 14)
#   Part D  one results PDF                   (was Step 15)

# ##########################################################################################
# PART A -- TEXT REPORT (was Step 12)
# ##########################################################################################
# =====================================================================
# REPORT: both cycles side by side (run after Steps 1-11)
# =====================================================================
REPORT.clear()
def _line(c='='): rep(c * 84)
def _row(name, value, unit=''): rep(f"  {name:<40}{value:>20} {unit}")
def _row2(name, vals, unit=''):
    rep(f"  {name:<40}" + "".join(f"{v:>18}" for v in vals) + f"  {unit}")
def _tot(d, key):                     # kWh -> MWh over a period
    return np.sum(d[key][:-1]) * dt / 1e3
def _each(fn, fmt):                   # apply fn to each cycle's results, format
    return [format(fn(R[ct], ct), fmt) for ct in cycles]

_line()
rep("  LUNAR POWER SYSTEM REPORT")
rep("  Solar tower -> regolith TES -> sCO2 Brayton cycle: simple vs recuperated")
rep("  Layout of Fig. 2: TES gas loop (sCO2) + power cycle (sCO2) + coolant loop (water)")
_line()

rep("\n[1] INPUTS (shared)")
_line('-')
_row("Lunar day / night", f"{t_day:.0f} / {t_night:.0f}", "h")
_row("Time step", f"{dt*60:.0f}", "min")
_row("Electrical demand / generator eff.", f"{W_elec_target:.1f} / {eta_gen:.2f}", "kWe")
_row("Parasitic: controls / heliostat drives", f"{W_aux_fixed:.1f} kWe / {w_helio_drive:.1f} W/m2")
_row("Gas loop: pressure / dp at design flow", f"{P_loop/1e6:.0f} MPa / {dp_loop*100:.1f} %")
_row("Gas loop: blower eff. / receiver rise", f"{eta_blower:.2f} / {dT_charge:.0f} K")
_row("Gas loop: TES HX gap / receiver wall gap", f"{dT_TES_HX:.0f} / {dT_wall:.0f}", "K")
_row("Receiver gas outlet, highest", f"{T_charge_out_max:.0f}", "K")
_row("Coolant loop: fluid / pressure", f"{COOLANT} / {P_cool_loop/1e6:.1f} MPa")
_row("Coolant: precooler gap / max water T", f"{dT_precool:.0f} / {T_cool_max:.0f}", "K")
_row("Coolant: loop dp / pump eff.", f"{dp_cool_loop/1e3:.0f} kPa / {eta_pump:.2f}")
_row("TES temperature range", f"{T_TES_max:.0f} -> {T_TES_min:.0f}", "K")
_row("Regolith type / density", f"{regolith_type} / {rho_TES:.0f}", "kg/m3")
_row("Regolith cp at T_max / T_min", f"{cp_regolith(T_TES_max):.0f} / {cp_regolith(T_TES_min):.0f}", "J/kg-K")
_row("TES mass margin", f"{TES_mass_margin:.3f}")
_row("TES insulation / ground temperature", f"{L_ins:.2f} m / {T_ground:.0f} K")
_row("Compressor inlet: lowest / highest", f"{T_comp_in:.0f} / {T_comp_in_max:.0f} K")
_row("Compressor inlet / outlet pressure", f"{P_low/1e6:.2f} / {P_high/1e6:.2f} MPa")
_row("Pressure drops: PHX / recup / cooler", f"{dp_PHX*100:.1f} / {dp_recup*100:.1f} / {dp_cooler*100:.1f}", "%")
_row("Compressor / turbine efficiency", f"{eta_comp:.2f} / {eta_turb:.2f}")
_row("Recuperator effectiveness", f"{eps_recup:.2f}")
_row("PHX approach (gas -> cycle sCO2)", f"{dT_PHX:.0f}", "K")
_row("Solar: G / optics (rho x dust x intercept)", f"{G_sun:.0f} W/m2 / {eta_opt:.3f}")
_row("Heliostat cosine factor: noon / sunrise", f"{float(cosine_factor(t_day/2)):.3f} / {float(cosine_factor(1e-6)):.3f}")
_row("Receiver: absorptance / emissivity", f"{alpha_rec:.2f} / {eps_rec:.2f}")
_row("Receiver: concentration", f"{C_ratio:.0f}x")
_row("Collector mode / latitude", f"{collector_mode} / {latitude_deg:.0f} deg")
_row("Habitat: crew / volume / wall area", f"{n_crew} / {V_hab:.0f} m3 / {A_hab:.0f} m2")
_row("Habitat: indoor T / MLI emittance", f"{T_hab:.0f} K / {eps_MLI:.3f}")
_row("Habitat: share of 100 kWe used inside", f"{f_elec_inside:.2f}")
_row("Habitat radiator / heat-pump hot side", f"{T_hab_rad:.0f} / {T_hp_cond:.0f} K")
_row("Radiator: emissivity / area margin", f"{emissivity:.2f} / {rad_area_margin:.2f}")
_row("Radiator view to surface / shield", f"{F_rad_surf:.2f} / {'yes' if shield else 'no'}")
_row("Surface temperature night / noon", f"{T_surf_night:.0f} / {T_surf_noon:.0f}", "K")

rep(f"\n[2] CYCLE STATES AT DESIGN POINT (TES fully charged, night, TIT {TIT_max:.0f} K, "
    f"compressor inlet {T_comp_in:.0f} K)")
_line('-')
_designs = {}
for ct in cycles:
    _c, _states = cycle_state_table(TIT_max, ct)
    _designs[ct] = _c
    rep(f"  {ct.capitalize()} cycle")
    rep(f"    {'State':<26}{'T (K)':>10}{'P (MPa)':>10}")
    for name, T, P in _states:
        rep(f"    {name:<26}{T:>10.1f}{P/1e6:>10.2f}")

rep("\n[3] DESIGN-POINT PERFORMANCE")
_line('-')
_row2("", [ct.capitalize() for ct in cycles])
_row2("Thermal efficiency", [f"{_designs[ct]['eta']*100:.2f}" for ct in cycles], "%")
_row2("Carnot limit", [f"{(1 - T_comp_in/TIT_max)*100:.1f}"] * len(cycles), "%")
_row2("Turbine pressure ratio", [f"{_c['P'][2 if ct == 'simple' else 3] / _c['P'][3 if ct == 'simple' else 4]:.2f}"
      for ct, _c in _designs.items()])
_row2("Back-work ratio (W_c / W_t)", [f"{_designs[ct]['w_c']/_designs[ct]['w_t']:.3f}" for ct in cycles])
_row2("Recuperator pinch (smallest hot-cold gap)",
      [f"{_designs[ct]['dT_pinch']:.1f}" if ct == 'recuperated' else '-' for ct in cycles], "K")
_row2("Electricity to habitat", _each(lambda r, ct: r['p_design']['W_elec'], '.1f'), "kWe")
_row2("Parasitic loads (night)", _each(lambda r, ct: r['p_design']['W_aux'], '.1f'), "kWe")
_row2("Gross generator output", _each(lambda r, ct: r['p_design']['W_gross'], '.1f'), "kWe")
_row2("Cycle shaft power (net)", _each(lambda r, ct: r['p_design']['W_net'], '.1f'), "kW")
_row2("CO2 mass flow", _each(lambda r, ct: r['p_design']['m_dot'], '.2f'), "kg/s")
_row2("Turbine power", _each(lambda r, ct: r['p_design']['W_t'], '.1f'), "kW")
_row2("Compressor power", _each(lambda r, ct: r['p_design']['W_c'], '.1f'), "kW")
_row2("Recuperator duty", _each(lambda r, ct: r['p_design']['Q_recup'], '.1f'), "kW")
_row2("Heat from TES", _each(lambda r, ct: r['p_design']['Q_in'], '.1f'), "kW")
_row2("Cooler heat rejection", _each(lambda r, ct: r['p_design']['Q_out'], '.1f'), "kW")
_row2("Precooler inlet temperature", _each(lambda r, ct: r['p_design']['T_cooler_in'], '.1f'), "K")
_row2("TES gas flow", _each(lambda r, ct: r['p_design']['m_gas'], '.2f'), "kg/s")
_row2("PHX gas in -> out", _each(lambda r, ct: f"{r['p_design']['T_gas_hot']:.0f} -> "
      f"{r['p_design']['T_gas_cold']:.0f}", 's'), "K")
_row2("Blower power", _each(lambda r, ct: r['p_design']['W_blower'], '.2f'), "kWe")
_row2("Coolant water flow", _each(lambda r, ct: r['p_design']['m_water'], '.2f'), "kg/s")
_row2("Precooler water in -> out", _each(lambda r, ct: f"{r['p_design']['T_w_cold']:.0f} -> "
      f"{r['p_design']['T_w_hot']:.0f}", 's'), "K")
_row2("Coolant pump power", _each(lambda r, ct: r['p_design']['W_pump'], '.2f'), "kWe")

rep("\n[4] SIZING")
_line('-')
_row2("", [ct.capitalize() for ct in cycles])
_row2("Paper method m_TES (Eq. 21, 24)", _each(lambda r, ct: r['m_TES_design']/1e3, '.1f'), "t")
_row2("Night-only minimum m_TES", _each(lambda r, ct: r['m_TES_night']/1e3, '.1f'), "t")
_row2("Design m_TES (with margin)", _each(lambda r, ct: r['m_TES']/1e3, '.1f'), "t")
_row2("TES volume (Eq. 25)", _each(lambda r, ct: r['V_TES'], '.0f'), "m3")
_row2("TES surface (cube)", _each(lambda r, ct: 6 * r['V_TES'] ** (2 / 3), '.0f'), "m2")
_row2("Collector area", _each(lambda r, ct: r['A_col'], '.0f'), "m2")
_row2("Radiator area for the TES floor (Step 6)", [f"{A_RAD_conservative[ct]:.0f}" for ct in cycles], "m2")
_row2("Radiator panel area installed (peak use)", [f"{A_RAD[ct]:.0f}" for ct in cycles], "m2")
_row2("  full power at night down to TES", [f"{CYC[ct]['T_TES_full_power']:.0f}" for ct in cycles], "K")
_row2(f"  {rad_sides}-sided panels, mass at {rad_areal_mass:.0f} kg/m2",
      [f"{A_RAD[ct] * rad_areal_mass / 1e3:.1f}" for ct in cycles], "t")
_row2("PHX size UA", [f"{UA_PHX[ct] / 1e3:.1f}" for ct in cycles], "kW/K")
_row2("Usable TES lower limit: cycle / plant",
      [f"{CYC[ct]['T_TES_cutoff_cycle']:.0f} / {CYC[ct]['T_TES_cutoff']:.0f}" for ct in cycles], "K")

rep("\n[5] OPERATION THROUGH THE LUNAR CYCLE")
_line('-')
_row2("", [ct.capitalize() for ct in cycles])
_row2("TES at sunrise", _each(lambda r, ct: r['day']['T_TES'][0], '.1f'), "K")
_row2("TES lowest", _each(lambda r, ct: r['cycle']['T_TES'].min(), '.1f'), "K")
_row2("  at (time after sunrise)",
      _each(lambda r, ct: r['cycle']['t'][np.argmin(r['cycle']['T_TES'])], '.0f'), "h")
_row2("TES at sunset", _each(lambda r, ct: r['night']['T_TES'][0], '.1f'), "K")
_row2("Cycle efficiency, max (running)", _each(lambda r, ct: np.nanmax(r['cycle']['eta'])*100, '.1f'), "%")
_row2("Cycle efficiency, min (running)", _each(lambda r, ct: np.nanmin(r['cycle']['eta'])*100, '.1f'), "%")
_row2("CO2 mass flow, max", _each(lambda r, ct: r['cycle']['m_dot'].max(), '.2f'), "kg/s")
_row2("Turbine inlet, min-max", _each(lambda r, ct: f"{np.nanmin(r['cycle']['TIT']):.0f}-"
      f"{np.nanmax(r['cycle']['TIT']):.0f}", 's'), "K")
_row2("TES - TIT (TES HX + PHX), min-max", _each(lambda r, ct: f"{np.nanmin(r['cycle']['dT_PHX_t']):.1f}-"
      f"{np.nanmax(r['cycle']['dT_PHX_t']):.1f}", 's'), "K")
_row2("Compressor inlet, night (min-max)",
      _each(lambda r, ct: f"{np.nanmin(r['night']['T_comp_in']):.0f}-{np.nanmax(r['night']['T_comp_in']):.0f}", 's'), "K")
_row2("Compressor inlet, day (min-max)",
      _each(lambda r, ct: (f"{np.nanmin(r['day']['T_comp_in']):.0f}-{np.nanmax(r['day']['T_comp_in']):.0f}"
                           if r['day']['running'].any() else 'off'), 's'), "K")
_row2("Hours compressor inlet above lowest",
      _each(lambda r, ct: np.sum(r['cycle']['running'][:-1] &
                                 (r['cycle']['T_comp_in'][:-1] > T_comp_in + 0.05)) * dt, '.0f'), "h")
_row2("TES gas flow, min-max", _each(lambda r, ct: f"{np.nanmin(r['cycle']['m_gas']):.2f}-"
      f"{np.nanmax(r['cycle']['m_gas']):.2f}", 's'), "kg/s")
_row2("Gas-loop dp, max (design = dp_loop)", _each(lambda r, ct: np.nanmax(r['cycle']['dp_loop_t']) * 100, '.2f'), "%")
_row2("Receiver temperature, min-max (day)", _each(lambda r, ct: f"{np.nanmin(r['day']['T_rec']):.0f}-"
      f"{np.nanmax(r['day']['T_rec']):.0f}", 's'), "K")
_row2("Receiver gas outlet, max", _each(lambda r, ct: np.nanmax(r['cycle']['T_rec_out']), '.0f'), "K")
_row2("Water to radiator, max", _each(lambda r, ct: np.nanmax(r['cycle']['T_w_hot']), '.0f'), "K")
_row2("Radiator use, max (running)", _each(lambda r, ct: np.nanmax(np.where(r['cycle']['running'],
      r['cycle']['rad_use'], np.nan)) * 100, '.0f'), "%")
_row2("Hours at reduced power (radiator limit)",
      _each(lambda r, ct: np.sum(r['cycle']['part_load'][:-1]) * dt, '.0f'), "h")
_row2("  lowest power in those hours",
      _each(lambda r, ct: (f"{np.min(r['cycle']['W_elec'][:-1][r['cycle']['part_load'][:-1]]):.0f}"
                           if r['cycle']['part_load'][:-1].any() else '-'), 's'), "kWe")
_row2("Hours cycle off", _each(lambda r, ct: np.sum(~r['cycle']['running'][:-1]) * dt, '.0f'), "h")
_row2("Demand served (energy)", _each(lambda r, ct: _tot(r['cycle'], 'W_elec') /
      (W_elec_target * t_lunar / 1e3) * 100, '.1f'), "%")

rep("\n[6] ENERGY TOTALS PER LUNAR CYCLE (MWh)")
_line('-')
_row2("", [ct.capitalize() for ct in cycles])
for key, name in [('Q_sun_in', 'Concentrated light on receiver'), ('Q_solar', 'Solar heat into gas loop (net)'), ('Q_curtail', 'Solar curtailed (TES full)'),
                  ('Q_in', 'Heat drawn by cycle'), ('Q_in_day', '  by day'), ('Q_in_night', '  by night'),
                  ('W_gross', 'Gross electricity generated'), ('W_aux', '  parasitic loads'),
                  ('W_elec', '  electricity to habitat'), ('W_unserved', 'Electricity NOT served'),
                  ('Q_out', 'Waste heat'),
                  ('Q_hab', '  to habitat'), ('Q_rad', '  to radiator')]:
    if key.endswith('_day'):
        vals = _each(lambda r, ct: _tot(r['day'], 'Q_in'), '.2f')
    elif key.endswith('_night'):
        vals = _each(lambda r, ct: _tot(r['night'], 'Q_in'), '.2f')
    else:
        vals = _each(lambda r, ct, k=key: _tot(r['cycle'], k), '.2f')
    _row2(name, vals, "MWh")
_row2("TES heat loss (insulation)", _each(lambda r, ct: _tot(r['cycle'], 'Q_loss'), '.2f'), "MWh")
_row2("  as % of solar heat to TES", _each(lambda r, ct: _tot(r['cycle'], 'Q_loss') /
      _tot(r['cycle'], 'Q_solar') * 100, '.1f'), "%")
_row2("  average / max", _each(lambda r, ct: f"{np.mean(r['cycle']['Q_loss']):.1f} / "
      f"{np.max(r['cycle']['Q_loss']):.1f}", 's'), "kW")
_row2("Receiver loss (reflection + radiation)", _each(lambda r, ct: _tot(r['cycle'], 'Q_rec_loss'), '.2f'), "MWh")
_row2("Receiver efficiency, lunar-day average",
      _each(lambda r, ct: _tot(r['cycle'], 'Q_solar') / _tot(r['cycle'], 'Q_sun_in') * 100, '.1f'), "%")
_row2("Receiver efficiency, at noon",
      _each(lambda r, ct: r['cycle']['eta_rec_t'][n_day // 2] * 100, '.1f'), "%")
_row2("Daytime hours receiver closed (loss > input)",
      _each(lambda r, ct: np.sum((r['day']['Q_sun_in'][:-1] > 0) & (r['day']['Q_solar'][:-1] <= 0)) * dt, '.0f'), "h")
_row2("  controls (always)", _each(lambda r, ct: W_aux_fixed * t_lunar / 1e3, '.2f'), "MWh")
_row2("  heliostat drives (day)", _each(lambda r, ct: _tot(r['cycle'], 'W_helio'), '.2f'), "MWh")
_row2("  gas-loop blower (day + night)", _each(lambda r, ct: _tot(r['cycle'], 'W_blower'), '.2f'), "MWh")
_row2("    of which by night", _each(lambda r, ct: _tot(r['night'], 'W_blower'), '.2f'), "MWh")
_row2("  coolant pump (day + night)", _each(lambda r, ct: _tot(r['cycle'], 'W_pump'), '.2f'), "MWh")
_row2("  parasitic load, night (min-max)", _each(lambda r, ct: f"{np.min(r['night']['W_aux']):.1f}-"
      f"{np.max(r['night']['W_aux']):.1f}", 's'), "kWe")
_row2("  peak parasitic load", _each(lambda r, ct: np.max(r['cycle']['W_aux']), '.1f'), "kWe")
_row2("Parasitic share of gross output",
      _each(lambda r, ct: _tot(r['cycle'], 'W_aux') / _tot(r['cycle'], 'W_gross') * 100, '.1f'), "%")
_row2("Light on receiver -> electricity",
      _each(lambda r, ct: _tot(r['cycle'], 'W_elec') / _tot(r['cycle'], 'Q_sun_in') * 100, '.2f'), "%")
_row2("Curtailed fraction of solar heat",
      _each(lambda r, ct: _tot(r['cycle'], 'Q_curtail') / _tot(r['cycle'], 'Q_solar') * 100, '.1f'), "%")

rep("\n[6b] HABITAT HEAT BALANCE (same for both cycles)")
_line('-')
_hb = R[cycles[0]]['cycle']
_row("Wall heat loss, night (max)", f"{np.max(_hb['Q_hab_wall']):.1f}", "kW")
_row("Wall heat gain, day (max)", f"{-np.min(_hb['Q_hab_wall']):.1f}", "kW")
_row("Heating needed from cycle, max", f"{np.max(_hb['Q_hab_need']):.1f}", "kW")
_row("Hours habitat needs heating", f"{np.sum(_hb['Q_hab_need'][:-1] > 0) * dt:.0f}", "h")
_row("Cooling load, min / max", f"{np.min(_hb['Q_hab_cool']):.1f} / {np.max(_hb['Q_hab_cool']):.1f}", "kW")
_row("Hours heat pump on", f"{np.sum(_hb['hp_on'][:-1]) * dt:.0f}", "h")
_row("Heat-pump power, max (inside 100 kWe)", f"{np.max(_hb['W_hp']):.1f}", "kWe")
_row("Habitat radiator panel area needed, max", f"{np.max(_hb['A_hab_rad']):.0f}", "m2")
_row("  mass at {:.0f} kg/m2".format(rad_areal_mass), f"{np.max(_hb['A_hab_rad']) * rad_areal_mass / 1e3:.1f}", "t")
_row2("", [ct.capitalize() for ct in cycles])
_row2("Heat sent to habitat from cycle", _each(lambda r, ct: _tot(r['cycle'], 'Q_hab'), '.2f'), "MWh")
_row2("Habitat heating not met", _each(lambda r, ct: _tot(r['cycle'], 'Q_hab_unmet'), '.2f'), "MWh")
_row2("Share of waste heat to habitat", _each(lambda r, ct: _tot(r['cycle'], 'Q_hab') /
      max(_tot(r['cycle'], 'Q_out'), 1e-9) * 100, '.1f'), "%")

rep("\n[7] CHECKS")
_line('-')
_row2("", [ct.capitalize() for ct in cycles])
_row2("Cycle energy balance residual", _each(lambda r, ct: r['cycle_balance'], '.1e'), "kW")
_row2("TES energy balance residual", _each(lambda r, ct: r['balance'], '.1e'), "kWh")
_row2("Cycle repeats (T_TES change)", _each(lambda r, ct: r['periodic_error'], '.3f'), "K")
for ct in cycles:
    r = R[ct]
    w = []
    if not r['survived']:
        w.append("TES ran out")
    if abs(r['periodic_error']) > 1.0:
        w.append("not yet repeating: raise n_cycles_periodic")
    _pl = np.sum(r['cycle']['part_load'][:-1]) * dt
    if _pl > 0 or (~r['cycle']['running'][:-1]).any():
        w.append(f"reduced power {_pl:.0f} h, {_tot(r['cycle'], 'W_unserved'):.1f} MWh not served: "
                 f"radiator cannot cope (raise rad_area_margin or T_comp_in_max, or set shield = True)")
    if r['cycle']['Q_hab_unmet'][:-1].sum() > 1e-6:
        w.append("habitat heating not fully met by the cycle waste heat")
    if np.nanmin(r['cycle']['eta']) < 0.08:
        w.append("cycle efficiency drops below 8%")
    if abs(r['cycle_balance']) > 1e-6 or abs(r['balance']) > 1e-3:
        w.append("energy balance residual not ~0")
    rep(f"  Status ({ct}): " + ("ALL OK" if not w else "WARNING - " + "; ".join(w)))
_line()

# ##########################################################################################
# PART B -- FIGURES (was Step 13)
# ##########################################################################################
# =====================================================================
# Step 13: Plots in AIAA journal format, simple vs recuperated (run after Steps 1-12)
# AIAA "Guidelines for Journal Figures and Tables":
#   - single-column width 3.25 in; lettering >= 8 pt at that width
#   - line weight >= 0.28 pt, distinct lines (here: colour AND line style, so the
#     figures also read in grayscale), no screens lighter than 20 %
#   - no title inside the figure: the caption goes centred BELOW the figure,
#     at most 20-25 words; axis labels in words with SI units ("Temperature, K")
#   - line art at >= 600 dpi
# Each figure is drawn at 3.25 x 2.6 in, Times-style serif font. Step 15 puts one
# figure per PDF page with its caption, and (SAVE_FIGS) also writes every figure
# at true size to FIG_DIR as vector PDF and 600-dpi PNG for the paper.
# figs = list of (name, figure, caption)
# =====================================================================
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

SAVE_FIGS = True                      # also write each figure to FIG_DIR (PDF + 600-dpi PNG)
FIG_DIR = 'figures_AIAA'
FIG_W, FIG_H = 3.25, 2.6              # in, AIAA single column

plt.rcdefaults()
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'Liberation Serif', 'Nimbus Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8, 'axes.titlesize': 8,
    'axes.linewidth': 0.6, 'lines.linewidth': 1.0, 'patch.linewidth': 0.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 3.5, 'ytick.major.size': 3.5, 'xtick.minor.size': 2, 'ytick.minor.size': 2,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'axes.grid': False, 'legend.frameon': False, 'legend.handlelength': 2.4,
    'figure.dpi': 150, 'savefig.dpi': 600, 'figure.max_open_warning': 0, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# --- One colour AND one line style per cycle (legible in colour and in grayscale) ---
COL = {'simple': '#1f4e9c', 'recuperated': '#b5281e'}
LS = {'simple': '-', 'recuperated': '--'}
LAB = {'simple': 'Simple', 'recuperated': 'Recuperated'}
REF = '0.35'                                  # reference lines: dark gray (>= 20 % screen)

days = R[cycles[0]]['cycle']['t'] / 24.0      # time since sunrise, Earth days
sunset_d = t_day / 24.0
XLAB_T = 'Time since sunrise, Earth days'

figs = []

def new_fig():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    return fig, ax

def time_axis(ax, day_only=False):
    """Time in Earth days, values every 5 days; sunset marked, 'Day'/'Night' labels."""
    end = sunset_d if day_only else days[-1]
    ax.set_xlim(0, end)
    ax.xaxis.set_major_locator(MultipleLocator(2 if day_only else 5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(XLAB_T)
    if not day_only:
        ax.axvline(sunset_d, color=REF, lw=0.6, ls=':')
        ax.text(sunset_d / 2, 1.02, 'Day', transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', fontsize=8)
        ax.text((sunset_d + days[-1]) / 2, 1.02, 'Night', transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', fontsize=8)

def both(ax, key, scale=1.0, **kw):
    for ct in cycles:
        ax.plot(days, R[ct]['cycle'][key] * scale, color=COL[ct], ls=LS[ct], label=LAB[ct], **kw)

def finish(fig, ax, name, ylabel, caption, legend_loc='best', ncol=1, headroom=0.0):
    ax.set_ylabel(ylabel)
    if headroom:                       # extra space above the data so the legend sits clear of the lines
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + (hi - lo) * headroom)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc=legend_loc, ncol=ncol)
    fig.tight_layout(pad=0.3)
    figs.append((name, fig, caption))

# 1. TES temperature -------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'T_TES')
ax.axhline(T_TES_max, color=REF, lw=0.6, ls=':')
for ct in cycles:
    ax.axhline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
ax.set_ylim(min(CYC[ct]['T_TES_full_power'] for ct in cycles) - 40, T_TES_max + 20)
time_axis(ax)
finish(fig, ax, 'TES_temperature', 'TES temperature, K',
       'TES temperature over one lunar cycle (dotted: maximum, and lowest TES giving full power with the installed radiator).',
       legend_loc='upper left', headroom=0.15)

# 2. Stored energy -----------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'E_TES', 1e-3)
time_axis(ax)
finish(fig, ax, 'TES_energy', 'Stored energy, MWh',
       f'Usable energy stored in the TES above {T_TES_min:.0f} K over one lunar cycle.', legend_loc='upper right')

# 3. Cycle efficiency ---------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'eta', 100)
time_axis(ax)
finish(fig, ax, 'cycle_efficiency', 'Cycle efficiency, %',
       'Cycle thermal efficiency (shaft power per unit heat from the TES) over one lunar cycle.',
       legend_loc='upper left', headroom=0.3)

# 4. Heat drawn from the TES -----------------------------------------------------------
fig, ax = new_fig()
both(ax, 'Q_in')
time_axis(ax)
finish(fig, ax, 'heat_drawn', 'Heat from TES, kW',
       f'Heat drawn from the TES to supply {W_elec_target:.0f} kWe to the habitat plus parasitic loads.',
       legend_loc='upper center', ncol=2, headroom=0.25)

# 5. Solar input vs demand -----------------------------------------------------------------
fig, ax = new_fig()
_ls2 = {'simple': ('-', '-.'), 'recuperated': ('--', ':')}
for ct in cycles:
    ax.plot(days, R[ct]['cycle']['Q_solar'], color=COL[ct], ls=_ls2[ct][0], label=f"{LAB[ct]}, solar collected")
    ax.plot(days, R[ct]['cycle']['Q_in'], color=COL[ct], ls=_ls2[ct][1], label=f"{LAB[ct]}, drawn by cycle")
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'solar_input', 'Thermal power, kW',
       f'Net solar heat into the gas loop and heat drawn by the cycle ({collector_mode} field, latitude {latitude_deg:.0f} deg).',
       legend_loc='upper right', headroom=0.5)

# 6. Parasitic loads, one figure per cycle, split by source ----------------------------------
_ymax = max(np.nanmax(R[ct]['cycle']['W_aux']) for ct in cycles) * 1.45
for ct in cycles:
    cy = R[ct]['cycle']
    fig, ax = new_fig()
    ax.plot(days, cy['W_aux'], color='k', lw=1.2, ls='-', label='Total')
    ax.plot(days, cy['W_blower'], color=COL[ct], ls='--', label='Gas-loop blower')
    ax.plot(days, cy['W_helio'], color='0.35', ls='-.', label='Heliostat drives')
    ax.plot(days, np.full_like(days, W_aux_fixed), color='0.35', ls=':', label='Controls')
    ax.plot(days, cy['W_pump'], color=COL[ct], ls=(0, (6, 2, 1, 2, 1, 2)), lw=0.9, label='Coolant pump')
    ax.set_ylim(0, _ymax)
    time_axis(ax)
    finish(fig, ax, f'parasitic_{ct}', 'Parasitic load, kWe',
           f'Parasitic electric load of the {ct} cycle plant by source; blower and pump run day and night.',
           legend_loc='upper right', ncol=2)

# 6b. Gas flow in the TES loop ---------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'm_gas')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'gas_flow', 'TES gas flow, kg/s',
       f'sCO2 flow in the TES gas loop: by day set by the heat stored (rise up to {dT_charge:.0f} K), at night by the PHX.',
       legend_loc='upper right', headroom=0.2)

# 6c. Gas-loop temperatures -------------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['T_gas_cold'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, PHX gas outlet')
    ax.plot(days, cy['T_rec_out'], color=COL[ct], ls='-.' if ct == 'simple' else ':',
            label=f'{LAB[ct]}, receiver gas outlet')
time_axis(ax)
finish(fig, ax, 'gas_loop_temperature', 'Temperature, K',
       'TES gas loop: gas leaving the PHX (back to receiver or TES) and leaving the receiver (day).',
       legend_loc='upper right', headroom=0.45)

# 6d. Coolant-loop temperatures ---------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['T_w_hot'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, from precooler')
    ax.plot(days, cy['T_w_cold'], color=COL[ct], ls='-.' if ct == 'simple' else ':',
            label=f'{LAB[ct]}, to precooler')
ax.axhline(T_cool_max, color=REF, lw=0.6, ls=':')
time_axis(ax)
finish(fig, ax, 'coolant_temperature', 'Water temperature, K',
       f'Coolant-loop water leaving and entering the precooler (dotted: {T_cool_max:.0f} K limit).',
       legend_loc='center right', headroom=0.1)

# 7. Compressor inlet ----------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'T_comp_in')
ax.axhline(T_comp_in_max, color=REF, lw=0.6, ls=':')
ax.set_ylim(T_comp_in - 5, T_comp_in_max + 5)
time_axis(ax)
finish(fig, ax, 'compressor_inlet', 'Compressor inlet temperature, K',
       f'Compressor inlet temperature set by the radiator (allowed {T_comp_in:.0f}-{T_comp_in_max:.0f} K, dotted: upper limit).',
       legend_loc='center right')

# 8. Waste heat to radiator ---------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'Q_rad')
time_axis(ax)
finish(fig, ax, 'waste_heat', 'Heat rejected, kW', 'Waste heat rejected by the plant radiator over one lunar cycle.',
       legend_loc='upper center', ncol=2, headroom=0.25)

# 9. Radiator use ---------------------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, np.where(cy['running'], cy['rad_use'], np.nan) * 100, color=COL[ct], ls=LS[ct], label=LAB[ct])
ax.axhline(100, color=REF, lw=0.6, ls=':')
ax.set_ylim(0, 105)
time_axis(ax)
finish(fig, ax, 'radiator_use', 'Radiator use, %',
       'Radiator area needed as a share of the installed area, sized to the peak hour ('
       + ', '.join(f"{LAB[ct].lower()} {A_RAD[ct]:,.0f} m$^2$" for ct in cycles) + ').',
       legend_loc='lower left')

# 10. Radiator environment -------------------------------------------------------------------
_cy = R[cycles[0]]['cycle']
fig, ax = new_fig()
ax.plot(days, T_surface(_cy['t']), color='k', lw=0.8, ls='-', label='Lunar surface')
ax.plot(days, _cy['T_sink'], color='0.45', lw=0.8, ls='-.', label='Radiator sink')
for ct in cycles:
    ax.plot(days, R[ct]['cycle']['T_comp_in'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, compressor inlet')
ax.set_ylim(0, max(T_surf_noon, T_comp_in_max) * 1.05)
time_axis(ax)
finish(fig, ax, 'radiator_environment', 'Temperature, K',
       'Lunar surface, radiator sink (' + ('shielded' if shield else f'view factor {F_rad_surf}')
       + ') and compressor inlet temperatures.', legend_loc='upper right', headroom=0.45)

# 11. TES and turbine inlet -----------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['TIT'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, turbine inlet')
    ax.plot(days, cy['T_TES'], color=COL[ct], ls='none', marker='o', ms=1.6, markevery=24,
            label=f'{LAB[ct]}, TES')
time_axis(ax)
finish(fig, ax, 'TES_and_TIT', 'Temperature, K',
       'Turbine inlet temperature (lines) and TES temperature (markers, every 24 h).', legend_loc='lower center')

# 12. PHX gap -------------------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'dT_PHX_t')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'PHX_gap', 'TES minus turbine inlet, K',
       f'TES minus turbine inlet: {dT_TES_HX + dT_PHX:.0f} K at night, {dT_PHX - dT_TES_HX:.0f} K while the sun charges the TES.',
       legend_loc='upper right', headroom=0.3)

# 13. Receiver efficiency (day) -------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, np.where(cy['Q_sun_in'] > 0, cy['eta_rec_t'], np.nan) * 100,
            color=COL[ct], ls=LS[ct], label=LAB[ct])
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_efficiency', 'Receiver efficiency, %',
       f'Receiver efficiency over the lunar day ({C_ratio:.0f}x concentration, emittance {eps_rec}); '
       f'receiver temperature from the gas loop.', legend_loc='best')

# 13b. Receiver temperature (day) -------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, np.where(cy['Q_sun_in'] > 0, cy['T_rec'], np.nan), color=COL[ct], ls=LS[ct], label=LAB[ct])
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_temperature', 'Receiver temperature, K',
       f'Receiver wall temperature: mean gas temperature in the receiver plus {dT_wall:.0f} K.', legend_loc='best')

# 14. Receiver loss (day) ------------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, np.where(cy['Q_sun_in'] > 0, cy['Q_rec_loss'], np.nan), color=COL[ct], ls=LS[ct], label=LAB[ct])
ax.set_ylim(bottom=0)
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_loss', 'Receiver loss, kW',
       'Receiver loss (reflection and thermal radiation) over the lunar day.', legend_loc='best', headroom=0.15)

# 15. Habitat heat balance ------------------------------------------------------------------
_hb = R[cycles[0]]['cycle']
fig, ax = new_fig()
ax.plot(days, _hb['Q_hab_cool'], color='k', ls='-', label='Cooling load')
ax.plot(days, _hb['Q_hab_wall'], color='0.45', ls='--', label='Wall loss (+) / gain (-)')
ax.plot(days, _hb['Q_hab_need'], color=COL['recuperated'], ls=':', label='Heating from cycle')
ax.plot(days, _hb['W_hp'], color=COL['simple'], ls='-.', label='Heat-pump power')
ax.axhline(0, color=REF, lw=0.4)
time_axis(ax)
finish(fig, ax, 'habitat_balance', 'Power, kW',
       f'Habitat heat balance ({n_crew} crew, {V_hab:.0f} m$^3$, {T_hab:.0f} K); the same for both cycles.',
       legend_loc='center right')

# 16. Habitat radiator area -----------------------------------------------------------------
fig, ax = new_fig()
ax.plot(days, _hb['A_hab_rad'], color='k', ls='-')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'habitat_radiator', 'Radiator area, m$^2$', 'Habitat radiator panel area needed over one lunar cycle.')

# 17-18. Performance versus TES temperature (night) ---------------------------------------
# Swept with the radiator NOT limiting (compressor inlet held at T_comp_in), so the curves
# show the cycles themselves over the whole TES window; dotted lines mark the lowest TES
# at which the installed radiator still gives full power.
_sweep, _A_keep = {}, dict(A_RAD)
for ct in cycles:
    A_RAD[ct] = np.inf
    T_sw = np.linspace(CYC[ct]['T_TES_cutoff'], T_TES_max, 160)
    P_sw = [plant_point(T, ct=ct) for T in T_sw]
    _sweep[ct] = (T_sw, [p if (p['running'] and not p['part_load']) else None for p in P_sw])
A_RAD.update(_A_keep)
def _sw(ct, key, scale=1.0):
    T_sw, P_sw = _sweep[ct]
    return T_sw, np.array([p[key] * scale if p is not None else np.nan for p in P_sw])
fig, ax = new_fig()
for ct in cycles:
    ax.plot(*_sw(ct, 'eta', 100), color=COL[ct], ls=LS[ct], label=LAB[ct])
    ax.axvline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
ax.set_xlabel('TES temperature, K')
ax.invert_xaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
finish(fig, ax, 'eta_vs_TES', 'Cycle efficiency, %',
       'Cycle efficiency versus TES temperature at night, radiator not limiting (dotted: full-power limit of installed radiator).',
       legend_loc='upper right')

fig, ax = new_fig()
for ct in cycles:
    ax.plot(*_sw(ct, 'Q_in'), color=COL[ct], ls=LS[ct], label=LAB[ct])
    ax.axvline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
ax.set_xlabel('TES temperature, K')
ax.invert_xaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
finish(fig, ax, 'Qin_vs_TES', 'Heat from TES, kW',
       f'Heat drawn from the TES for {W_elec_target:.0f} kWe at night, radiator not limiting (dotted: full-power limit of installed radiator).',
       legend_loc='upper left')

# 19. T-s diagram ---------------------------------------------------------------------------
def s_of(h, P):
    return CP.PropsSI('S', 'H', h, 'P', P, FLUID) / 1e3
def T_of(h, P):
    return CP.PropsSI('T', 'H', h, 'P', P, FLUID)
def isobar(P, ha, hb, n=60):
    h = np.linspace(ha, hb, n)
    return np.array([s_of(x, P) for x in h]), np.array([T_of(x, P) for x in h])

fig, ax = new_fig()
Ts = np.linspace(220, 304.0, 120)
ax.plot([CP.PropsSI('S', 'T', T, 'Q', 0, FLUID) / 1e3 for T in Ts], Ts, color='0.45', lw=0.6)
ax.plot([CP.PropsSI('S', 'T', T, 'Q', 1, FLUID) / 1e3 for T in Ts], Ts, color='0.45', lw=0.6,
        label='Saturation dome')
for ct in cycles:
    if ct == 'simple':
        r = simple_cycle(TIT_max)
        pts = list(zip([r['h1'], r['h2'], r['h3'], r['h4']], r['P']))
        legs = [('line', 0, 1), ('isobar', 1, 2), ('line', 2, 3), ('isobar', 3, 0)]
    else:
        r = recuperated_cycle(TIT_max)
        pts = list(zip([r['h1'], r['h2'], r['h3'], r['h4'], r['h5'], r['h6']], r['P']))
        legs = [('line', 0, 1), ('isobar', 1, 2), ('isobar', 2, 3), ('line', 3, 4),
                ('isobar', 4, 5), ('isobar', 5, 0)]
    for kind, i, j in legs:
        (hi, Pi), (hj, Pj) = pts[i], pts[j]
        if kind == 'line':
            ax.plot([s_of(hi, Pi), s_of(hj, Pj)], [T_of(hi, Pi), T_of(hj, Pj)], color=COL[ct], ls=LS[ct])
        else:
            x, y = isobar(Pi, hi, hj)
            ax.plot(x, y, color=COL[ct], ls=LS[ct])
    ax.plot([], [], color=COL[ct], ls=LS[ct], label=f"{LAB[ct]} ({r['eta']*100:.1f} %)")
    for n, (h, P) in enumerate(pts, 1):
        ax.plot(s_of(h, P), T_of(h, P), 'o', color=COL[ct], ms=2.5)
        ax.annotate(str(n), (s_of(h, P), T_of(h, P)), textcoords='offset points',
                    xytext=(3, -8) if ct == 'recuperated' else (-7, 3), color=COL[ct], fontsize=8)
ax.set_xlabel('Specific entropy, kJ/(kg$\\cdot$K)')
ax.xaxis.set_minor_locator(AutoMinorLocator())
finish(fig, ax, 'Ts_diagram', 'Temperature, K',
       f'Temperature-entropy diagram at full charge (turbine inlet {TIT_max:.0f} K); numbers mark the cycle states.',
       legend_loc='upper left')

# 20. Energy totals ------------------------------------------------------------------------
items = [('Q_sun_in', 'Light'), ('Q_solar', 'Col-\nlected'), ('Q_loss', 'TES\nloss'),
         ('Q_in', 'To\ncycle'), ('W_aux', 'Para-\nsitic'), ('W_elec', 'To\nhabitat'),
         ('Q_rad', 'To\nradiator')]
tot = lambda d, k: np.sum(d[k][:-1]) * dt / 1e3
x = np.arange(len(items)); w = 0.38
fig, ax = new_fig()
_hatch = {'simple': '', 'recuperated': '////'}
for i, ct in enumerate(cycles):
    off = (-1 if i == 0 else 1) * (w / 2 + 0.01)
    ax.bar(x + off, [tot(R[ct]['cycle'], k) for k, _ in items], w, color='white' if i else COL[ct],
           edgecolor=COL[ct], hatch=_hatch[ct], lw=0.6, label=LAB[ct])
ax.set_xticks(x, [lbl for _, lbl in items])
ax.tick_params(axis='x', which='both', top=False, bottom=False)
ax.xaxis.set_minor_locator(plt.NullLocator())
ax.set_ylim(bottom=0)
finish(fig, ax, 'energy_totals', 'Energy per lunar cycle, MWh',
       'Energy per lunar cycle: light on receiver, solar heat collected, TES loss, heat to cycle, parasitic, habitat, radiator.',
       legend_loc='upper right')

# --- Nothing shown in the cell; figures go into the PDF (Part D) and, if SAVE_FIGS, to FIG_DIR ---
# Each figure is rendered to a 600-dpi PNG ONCE (FIG_PNG); the same image is written to
# FIG_DIR and placed on its PDF page in Part D (before, it was rendered twice).
import os
FIG_PNG = []
if SAVE_FIGS:
    os.makedirs(FIG_DIR, exist_ok=True)
for _k, (_n, _f, _c) in enumerate(figs, 1):
    _buf = io.BytesIO()
    _f.savefig(_buf, format='png', dpi=600)
    FIG_PNG.append(_buf.getvalue())
    if SAVE_FIGS:
        _f.savefig(os.path.join(FIG_DIR, f'Fig{_k:02d}_{_n}.pdf'))
        with open(os.path.join(FIG_DIR, f'Fig{_k:02d}_{_n}.png'), 'wb') as _fh:
            _fh.write(FIG_PNG[-1])
for _n, _f, _c in figs:
    plt.close(_f)
say(f"Step 13: {len(figs)} AIAA-format figures ({FIG_W} x {FIG_H} in, 8 pt serif)"
    + (f", saved to {FIG_DIR}/ as PDF and 600-dpi PNG" if SAVE_FIGS else ''))

# ##########################################################################################
# PART C -- REFERENCES (was Step 14)
# ##########################################################################################
# =====================================================================
# REFERENCES FOR THE INPUT VALUES (Step 2), for citing in the paper
# Each input is tagged in Step 2 with its reference number [R..].
# 'Assumed' = engineering assumption, no single source: state it in the paper
# and cover it with a sensitivity study.
# =====================================================================
REFERENCES = [
    # tag,  inputs,  value used,  source
    ("R1", "G_sun", "1361 W/m2",
     "Kopp, G., Lean, J.L. (2011). A new, lower value of total solar irradiance: Evidence and "
     "climate significance. Geophysical Research Letters 38, L01706. doi:10.1029/2010GL045777"),
    ("R2", "rho_mirror", "0.95",
     "Typical silvered-glass heliostat reflectance (~95 %); assumed, see e.g. Stine, W.B., "
     "Geyer, M. (2001). Power From The Sun, Ch. 10 (central receivers)"),
    ("R3", "f_dust", "0.93",
     "Terrestrial heliostat soiling losses of 6.3-8.8 % (Strachan et al., 1993, Sandia), as "
     "summarised in Guan, Z. et al., Heat Exchanger Fouling and Cleaning conf. Lunar dust "
     "loss is not known: assumed"),
    ("R4", "alpha_rec, eps_rec", "0.95, 0.85",
     "Ho, C.K., Iverson, B.D. (2014). Review of high-temperature central receiver designs for "
     "concentrating solar power. Renewable and Sustainable Energy Reviews 29, 835-846. "
     "doi:10.1016/j.rser.2013.08.099"),
    ("R5", "C_ratio", "500 (towers 500-1000)",
     "Ho & Iverson (2014) [R4]: baseline concentration ~900; central receivers typically "
     "500-1000 suns. Ho, C.K. (2017). Advances in central receivers for concentrating solar "
     "applications. Solar Energy 152, 38-56. doi:10.1016/j.solener.2017.03.048"),
    ("R6", "cosine factor (Step 8)", "cos(0.5 * angle sun-receiver)",
     "Heliostat normal bisects the sun and receiver directions; effective area reduced by the "
     "cosine of half that angle: Stine & Geyer (2001) [R2], Ch. 10"),
    ("R7", "Regolith cp(T), k(T), rho_TES", "Table 3 correlations, 1500 kg/m3",
     "Paper Table 3, ref. [16]: Schreiner, S.S., Dominguez, J.A., Sibille, L., Hoffman, J.A. "
     "(2016). Thermophysical property models for lunar regolith. Advances in Space Research "
     "57(5), 1209-1222 -- check volume/pages/DOI against the paper's ref. [16]"),
    ("R8", "T_surf_noon, T_surf_night", "392 K, 100 K",
     "Vasavada, A.R. et al. (2012). Lunar equatorial surface temperatures and regolith "
     "properties from the Diviner Lunar Radiometer Experiment. JGR Planets 117, E00H18. "
     "doi:10.1029/2011JE003987 (385 K noon, 101 K midnight, 95 K pre-dawn); Williams, J.-P. "
     "et al. (2017). The global surface temperatures of the Moon as measured by Diviner. "
     "Icarus 283, 300-325. doi:10.1016/j.icarus.2016.08.012 (387-397 K at the equator)"),
    ("R9", "T_ground", "250 K",
     "Hayne, P.O. et al. (2017). Global regolith thermophysical properties of the Moon from "
     "Diviner. JGR Planets 122, 2371-2400. doi:10.1002/2017JE005387 (diurnal mean ~252-256 K "
     "at ~1 m depth, from Keihm et al. 1973 Apollo heat-flow data)"),
    ("R10", "eta_comp, eta_turb", "0.85, 0.90",
     "NETL (2014). An Assessment of Supercritical CO2 Power Cycles Integrated with Generic Heat "
     "Sources, DOE/NETL-4001/09/091014 (compressor 0.85); turbine 0.90 is a common sCO2 "
     "assumption (e.g. TESPy sCO2 recompression model, tespy.readthedocs.io)"),
    ("R11", "eps_recup", "0.95",
     "Effect of actual recuperators' effectiveness on the attainable efficiency of sCO2 Brayton "
     "cycles for solar thermal power plants. Energies 15 (2022) 7773. doi:10.3390/en15207773 "
     "(effectiveness 0.95 assumed for recuperators)"),
    ("R12", "rad_sides, rad_areal_mass, T_hab_rad", "2-sided, 14 kg/m2, 283 K",
     "Schunk, R.G. et al. (2022). Thermal Control System Architecture and Technology Challenges "
     "for a Lunar Surface Habitat. IEEE Aerospace Conf., NTRS 20210026557 (48 m2 double-sided "
     "radiators, ISS HRS areal density 14 kg/m2); Schunk et al., TFAWS 2021, NTRS 20210020015 "
     "(radiator temperature 283 K)"),
    ("R13", "V_hab, n_crew", "419 m3, 4 crew",
     "Kennedy, K.J. et al. (2008). Lunar Architecture Team Phase 2 Habitat Volume Estimation: "
     "'Caution When Using Analogs'. Earth & Space 2008, ASCE. NTRS 20080013517 (mean estimated "
     "total pressurized volume 419 m3 for 4 crew)"),
    ("R14", "eps_MLI", "0.02",
     "MLI effective emissivity typically 0.015-0.03: Applied Sciences 12 (2022) 2858. "
     "doi:10.3390/app12062858"),
    ("R15", "f_carnot_hp, T_hp_cond", "0.43, 360 K",
     "Goldman, J., Lovell, T., Ewert, M., Petete, P. (1994). High Temperature Lift Heat Pump "
     "Refrigerant and Thermodynamic Cycle Selection. SAE 941272 (COP 43 % of Carnot; optimum "
     "rejection 377-387 K); Sridhar & Gottmann, NASA CR-190063 (1992) (radiator ~362 K)"),
    ("R16", "q_metab", "14 MJ/person/day",
     "A More Completely Defined CELSS. SAE Technical Paper 941292 (1994) (metabolic load "
     "14 MJ/person/day)"),
    ("R17", "eta_gen, T_comp_in range, P ranges", "0.95, 305-320 K, 7.5-8.0 / 20-25 MPa",
     "Paper Table 2"),
    ("R19", "w_helio_drive", "0.5 W/m2",
     "Heliostat drive power: commercial fields have used ~0.5-4.5 % of plant output for "
     "tracking drives (US Patent 9,170,032, heliostat tracking and operation); 0.5 W/m2 is an "
     "assumed low-power value for small lunar heliostats"),
    ("R18", "PHX / gas-to-sCO2 heat transfer (context for dT_PHX)", "U ~ 300-400 W/m2-K",
     "Albrecht, K.J., Ho, C.K. et al. (Sandia), Development and testing of a 20 kWt moving "
     "packed-bed particle-to-sCO2 heat exchanger, OSTI 1883488"),
]

ASSUMED = [
    ("Layout of paper Fig. 2", "TES gas loop (sCO2) + power cycle + water coolant loop; receiver bypassed at night"),
    ("P_loop = 25 MPa, dp_loop = 3 % at the design flow", "TES gas loop; dp scales with m_gas^2/rho (design = largest of operation)"),
    ("eta_blower = 0.7", "TES gas-loop blower (runs day and night)"),
    ("dT_charge = 150 K", "largest gas rise above the TES gas exit while charging (sets the daytime gas flow)"),
    ("dT_TES_HX = 10 K, dT_wall = 20 K", "TES gas exit above/below T_TES; receiver wall above the mean gas temperature"),
    ("dT_PHX = 10 K", "PHX (gas -> cycle sCO2): hot-end approach and smallest gap"),
    ("Coolant: water, 1.5 MPa, <= 450 K", "coolant loop fluid and limits (may change)"),
    ("dT_precool = 10 K, dp_cool_loop = 150 kPa, eta_pump = 0.7", "precooler gap and coolant pump"),
    ("f_intercept = 0.95, beta_target_deg = 30", "heliostat spillage/blocking, field geometry"),
    ("dp_PHX, dp_recup, dp_cooler = 1 %", "component pressure drops (cycle side)"),
    ("Recuperator always in the flow path", "heat flows backwards when the turbine exhaust is colder than the compressor outlet"),
    ("W_aux_fixed = 3 kWe", "controls, sensors, power electronics"),
    ("TES_mass_margin = 1.20, rad_area_margin = 1.00", "design margins; radiator = peak need of normal operation"),
    ("L_ins = 1 m", "loose-regolith insulation around the TES"),
    ("D_hab = 6 m, f_elec_inside = 1.0", "habitat shape; all 100 kWe dissipated inside"),
    ("shield = True, F_rad_surf = 0.5", "radiators shaded from the lunar surface"),
    ("radiator fin efficiency = 1", "panel at the local water temperature (real panels 0.8-0.9)"),
    ("f_intercept constant", "no extra heliostat shading/blocking at low sun (optimistic at high latitude)"),
    ("T_comp_in_max = 380 K", "highest compressor inlet before reduced power"),
]

# Validation of the cycle model (Steps 4-5) against measured sCO2 Brayton data
VALIDATION = [
    ("V1", "Wright, S.A., Conboy, T.M., Rochau, G.E. (2012). High-Temperature Split-Flow Recompression "
     "Brayton Cycle Initial Test Results. SAND2012-6426, Sandia National Laboratories. OSTI 1051730. "
     "Run GenIV_110714_0952 (simple recuperated loop): state table at 5000 s used for validation"),
    ("V2", "Pasch, J., Carlson, M., Fleming, D., Rochau, G. (2016). Evaluation of Recent Data from the "
     "Sandia National Laboratories Closed Brayton Cycle Testing. ASME Turbo Expo GT2016-57620"),
    ("V3", "Conboy, T., Wright, S., Pasch, J., Fleming, D., Rochau, G., Fuller, R. (2012). Performance "
     "Characteristics of an Operating Supercritical CO2 Brayton Cycle. J. Eng. Gas Turbines Power "
     "134(11), 111703. doi:10.1115/1.4007199"),
    ("V4", "Wright, S.A., Radel, R.F., Vernon, M.E., Rochau, G.E., Pickard, P.S. (2010). Operation and "
     "Analysis of a Supercritical CO2 Brayton Cycle. SAND2010-0171 (near-critical compression)"),
    ("V5", "SAND2016-6214B: Sandia recompression loop design values (turbines ~85 %, main compressor "
     "66.5 %, recompressor 70.1 %)"),
]
VALIDATION_SUMMARY = [
    "Sandia run GenIV_110714_0952 (SAND2012-6426), rebuilt with the Step 4/5 equations at the measured",
    "boundary conditions (306 K / 589 K, 2.61 kg/s, 7.64 -> 12.24 MPa):",
    "  with Sandia design efficiencies: state temperatures within ~1 %, turbine power +3.7 %,",
    "  recuperator duty -2.1 %, heater duty +4.9 %; compressor power -18 % (as-tested efficiency",
    "  0.576 vs design 0.701 at a near-critical, part-speed point).",
    "  with the as-tested efficiencies: turbine, compressor and net shaft power reproduced exactly,",
    "  recuperator and heater duty within 1 %.",
    "Verification: energy balance closes to machine precision, entropy generated >= 0, below Carnot,",
    "helium ideal-gas limit matches the closed-form Brayton cycle within 0.04 %.",
]

REFS_LINES = []
def _ref_out(*args):
    REFS_LINES.extend(" ".join(str(a) for a in args).split("\n"))

_ref_out("REFERENCES FOR INPUT VALUES")
_ref_out("=" * 100)
for tag, inputs, value, source in REFERENCES:
    _ref_out(f"[{tag}] {inputs} = {value}")
    _ref_out(f"      {source}\n")
_ref_out("ASSUMED VALUES (state in the paper, cover with sensitivity)")
_ref_out("-" * 100)
for inp, why in ASSUMED:
    _ref_out(f"  {inp:<58} {why}")

_ref_out("")
_ref_out("VALIDATION OF THE CYCLE MODEL")
_ref_out("-" * 100)
for line in VALIDATION_SUMMARY:
    _ref_out(f"  {line}")
_ref_out("")
for tag, source in VALIDATION:
    _ref_out(f"[{tag}] {source}\n")

say(f"Step 14: {len(REFERENCES)} references, {len(ASSUMED)} assumptions and "
    f"{len(VALIDATION)} validation sources listed for the PDF")

# ##########################################################################################
# PART D -- RESULTS PDF (was Step 15)
# ##########################################################################################
# =====================================================================
# ONE RESULTS PDF: summary, run log, full report, all figures, references
# (run after Steps 1-14; nothing else is saved and nothing is printed in the cells)
# =====================================================================
import datetime

def _text_pages(pdf, title, lines, per_page=80, width=125, fontsize=7.0):
    """Write a list of text lines as monospace pages. Returns the number of pages."""
    wrapped = []
    for line in lines:
        pieces = textwrap.wrap(line, width, subsequent_indent='        ',
                               drop_whitespace=False, replace_whitespace=False)
        wrapped.extend(pieces if pieces else [''])
    n = 0
    for k in range(0, max(len(wrapped), 1), per_page):
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.06, 0.965, title + ('  (continued)' if k else ''), fontsize=11,
                 weight='bold')
        fig.text(0.06, 0.945, '\n'.join(wrapped[k:k + per_page]), fontsize=fontsize,
                 family='monospace', va='top')
        pdf.savefig(fig)
        plt.close(fig)
        n += 1
    return n

def _summary_lines():
    s = [f"Run: {datetime.datetime.now():%Y-%m-%d %H:%M}",
         f"Demand {W_elec_target:.0f} kWe | TES {T_TES_max:.0f} -> {T_TES_min:.0f} K | "
         f"compressor inlet {T_comp_in:.0f}-{T_comp_in_max:.0f} K at {P_low/1e6:.1f} MPa | "
         f"P_high {P_high/1e6:.0f} MPa | latitude {latitude_deg:.0f} deg",
         "", f"{'':42}" + "".join(f"{ct.capitalize():>18}" for ct in cycles)]
    def row(name, fn, fmt, unit=''):
        s.append(f"{name:<42}" + "".join(f"{format(fn(ct), fmt):>18}" for ct in cycles) + f"  {unit}")
    row("Cycle efficiency at full charge", lambda ct: cycle_at(TIT_max, ct)['eta'] * 100, '.1f', '%')
    row("TES mass", lambda ct: R[ct]['m_TES'] / 1e3, ',.0f', 't')
    row("TES volume", lambda ct: R[ct]['V_TES'], ',.0f', 'm3')
    row("Collector (mirror) area", lambda ct: R[ct]['A_col'], ',.0f', 'm2')
    row("Plant radiator panel area", lambda ct: A_RAD[ct], ',.0f', 'm2')
    row("Plant radiator mass", lambda ct: A_RAD[ct] * rad_areal_mass / 1e3, '.1f', 't')
    row("PHX size UA", lambda ct: UA_PHX[ct] / 1e3, '.1f', 'kW/K')
    row("Blower power, night average", lambda ct: np.mean(R[ct]['night']['W_blower']), '.1f', 'kWe')
    row("Parasitic load, day / night (avg)", lambda ct: f"{np.mean(R[ct]['day']['W_aux']):.1f} / "
        f"{np.mean(R[ct]['night']['W_aux']):.1f}", 's', 'kWe')
    row("TES lowest in operation", lambda ct: R[ct]['cycle']['T_TES'].min(), '.0f', 'K')
    row("Parasitic share of gross output", lambda ct: np.sum(R[ct]['cycle']['W_aux'][:-1]) /
        np.sum(R[ct]['cycle']['W_gross'][:-1]) * 100, '.1f', '%')
    row("Demand served (energy)", lambda ct: np.sum(R[ct]['cycle']['W_elec'][:-1]) * dt /
        (W_elec_target * t_lunar) * 100, '.1f', '%')
    row("Light on receiver -> electricity", lambda ct: np.sum(R[ct]['cycle']['W_elec'][:-1]) /
        np.sum(R[ct]['cycle']['Q_sun_in'][:-1]) * 100, '.1f', '%')
    s += ["", "Status lines from the report:"]
    s += [l.strip() for l in REPORT if l.strip().startswith('Status')]
    s += ["", "Contents: summary, run log, full report, figures, references and validation."]
    return s

pages = 0
with PdfPages(OUT_PDF) as pdf:
    pages += _text_pages(pdf, 'Lunar solar-TES sCO2 power system: summary', _summary_lines(),
                         fontsize=8.0, width=110)
    pages += _text_pages(pdf, 'Run log (messages from Steps 2-14)', RUN_LOG)
    pages += _text_pages(pdf, 'Report (Step 12)', REPORT)
    for k, (name, f, caption) in enumerate(figs, 1):
        # one figure per page: the AIAA-size figure (3.25 in) shown at 2x, caption centred below
        img = plt.imread(io.BytesIO(FIG_PNG[k - 1]))          # rendered once in Part B
        page = plt.figure(figsize=(8.5, 11))
        w_in = 2 * FIG_W
        h_in = w_in * img.shape[0] / img.shape[1]
        ax_img = page.add_axes([(8.5 - w_in) / 2 / 8.5, 1 - (1.0 + h_in) / 11, w_in / 8.5, h_in / 11])
        ax_img.imshow(img, interpolation='lanczos')
        ax_img.axis('off')
        page.text(0.5, 1 - (1.25 + h_in) / 11, f"Fig. {k}  {caption}", ha='center', va='top',
                  fontsize=11, family='serif', weight='bold', wrap=True)
        pdf.savefig(page)
        plt.close(page)
        pages += 1
    pages += _text_pages(pdf, 'References, assumptions and validation (Step 14)', REFS_LINES)
    info = pdf.infodict()
    info['Title'] = 'Lunar solar-TES sCO2 power system results'
    info['Author'] = 'Lunar_Regolith_TES_Cycle notebook'

print(f"Saved {OUT_PDF}: {pages} pages (summary, run log, report, {len(figs)} figures, "
      f"references and validation)")
