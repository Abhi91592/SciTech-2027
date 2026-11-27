# lunar_report.py -- report, figures, references and the results PDF, in one go
# Loaded by Step 12 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses the results of Steps 1-11 (R, CYC, A_RAD, ...).
#   Part A  text report (short: the numbers for the paper)
#   Part B  AIAA figures, saved to FIG_DIR
#   Part C  references, assumptions, validation
#   Part D  one results PDF: summary page, report, figures, references
# Time in the report and the figures = time since sunrise (the end of the longest
# darkness of the design month, Step 8). "Dark" = sun hidden at the middle of the step.

# --- What goes into the PDF -------------------------------------------------------------
PDF_RUN_LOG = False          # True: also add the run log (every message from Steps 2-11)
PDF_REFERENCES = True        # references, assumptions and validation pages

# ##########################################################################################
# PART A -- TEXT REPORT
# ##########################################################################################
REPORT.clear()
def _line(c='='): rep(c * 84)
def _row(name, value, unit=''): rep(f"  {name:<44}{value:>32} {unit}")
def _row2(name, vals, unit=''):
    rep(f"  {name:<44}" + "".join(f"{v:>16}" for v in vals) + f"  {unit}")
def _head(title):
    rep(f"\n{title}")
    _line('-')
def _cols():
    _row2("", [ct.capitalize() for ct in cycles])
def _tot(d, key):                     # kWh -> MWh over a period (steps, not the end point)
    return np.sum(d[key][:-1]) * dt / 1e3
def _each(fn, fmt):                   # fn(results of a cycle, cycle) for each cycle, formatted
    return [format(fn(R[ct], ct), fmt) for ct in cycles]
def _rng(a, fmt='.0f'):               # "min-max" of an array, ignoring NaN
    return f"{np.nanmin(a):{fmt}}-{np.nanmax(a):{fmt}}" if np.isfinite(a).any() else '-'

_cy0 = R[cycles[0]]['cycle']
_lit0 = _cy0['is_day'][:-1]
_cos_lit = float(np.mean(cosine_factor(_cy0['t'][:-1][_lit0] + 0.5 * dt))) if _lit0.any() else float('nan')

_line()
rep("  LUNAR POWER SYSTEM REPORT")
rep("  Solar tower -> regolith TES -> sCO2 Brayton cycle (simple vs recuperated)")
rep(f"  Site: {SITE.replace('_', ' ')}, latitude {latitude_deg:.2f} deg | time step {dt:g} h")
_line()

_head("[1] SITE AND SUN")
_row("Sunlit over one draconic year", f"{100*SITE_STATS['illum_fraction']:.1f}", "%")
_row("Longest darkness over the year", f"{SITE_STATS['longest_dark_h']:.0f}", "h")
_row("Sun elevation, min / max", f"{SITE_STATS['el_min_deg']:.1f} / {SITE_STATS['el_max_deg']:.1f}", "deg")
_row("Design month (darkest of the year)", f"{t_lunar:.0f}", "h")
_row("  sunlit / longest darkness", f"{t_day:.0f} / {t_night:.0f}", "h")
_row("  dark in total (all dark spells)", f"{SITE_STATS['dark_hours_month']:.0f}", "h")
_row("Heliostat field / mean cosine factor (lit)", f"{field_layout} / {_cos_lit:.3f}")

_head("[2] MAIN INPUTS (all inputs: lunar_inputs.py)")
_row("Electrical demand / generator efficiency", f"{W_elec_target:.0f} kWe / {eta_gen:.2f}")
_row("TES range / regolith / mass margin", f"{T_TES_max:.0f}-{T_TES_min:.0f} K / {regolith_type} / {TES_mass_margin:.2f}")
_row("Compressor inlet range / pressures", f"{T_comp_in:.0f}-{T_comp_in_max:.0f} K / {P_low/1e6:.0f}-{P_high/1e6:.0f} MPa")
_row("Compressor / turbine / recuperator eff.", f"{eta_comp:.2f} / {eta_turb:.2f} / {eps_recup:.2f}")
_row("Receiver: concentration / emittance / optics", f"{C_ratio:.0f}x / {eps_rec:.2f} / {eta_opt:.3f}")
_row("Radiator: emissivity / faces / mass", f"{emissivity:.2f} / {rad_sides} / {rad_areal_mass:.0f} kg/m2")

rep(f"\n[3] DESIGN POINT (TES full, no sun, turbine inlet {TIT_max:.0f} K, compressor inlet {T_comp_in:.0f} K)")
_line('-')
_designs = {ct: cycle_state_table(TIT_max, ct) for ct in cycles}
_cols()
_row2("Cycle efficiency", [f"{_designs[ct][0]['eta']*100:.1f}" for ct in cycles], "%")
_row2("  Carnot limit", [f"{(1 - T_comp_in/TIT_max)*100:.1f}"] * len(cycles), "%")
_row2("Heat from TES", _each(lambda r, ct: r['p_design']['Q_in'], '.0f'), "kW")
_row2("Shaft power (habitat + parasitic)", _each(lambda r, ct: r['p_design']['W_net'], '.1f'), "kW")
_row2("  parasitic loads", _each(lambda r, ct: r['p_design']['W_aux'], '.1f'), "kWe")
_row2("CO2 flow", _each(lambda r, ct: r['p_design']['m_dot'], '.2f'), "kg/s")
_row2("Recuperator duty", _each(lambda r, ct: r['p_design']['Q_recup'], '.0f'), "kW")
_row2("Waste heat", _each(lambda r, ct: r['p_design']['Q_out'], '.0f'), "kW")
rep("  Cycle states (T, P):")
for ct in cycles:
    _st = [f"{name.split()[0]} {T:.0f} K {P/1e6:.1f} MPa" for name, T, P in _designs[ct][1]]
    for _k in range(0, len(_st), 3):
        rep(f"    {ct.capitalize() if _k == 0 else '':<13}" + " | ".join(_st[_k:_k + 3]))

_head("[4] SIZING")
_cols()
_row2("TES mass, paper method (Eq. 21, 24)", _each(lambda r, ct: r['m_TES_design']/1e3, ',.0f'), "t")
_row2("TES mass, minimum for the longest darkness", _each(lambda r, ct: r['m_TES_night']/1e3, ',.0f'), "t")
_row2("TES mass, design (with margin)", _each(lambda r, ct: r['m_TES']/1e3, ',.0f'), "t")
_row2("TES volume", _each(lambda r, ct: r['V_TES'], ',.0f'), "m3")
_row2("Collector (mirror) area", _each(lambda r, ct: r['A_col'], ',.0f'), "m2")
_row2("Radiator panel area (peak use)", [f"{A_RAD[ct]:,.0f}" for ct in cycles], "m2")
_row2("Radiator mass", [f"{A_RAD[ct] * rad_areal_mass / 1e3:.1f}" for ct in cycles], "t")
_row2("PHX size UA", [f"{UA_PHX[ct] / 1e3:.1f}" for ct in cycles], "kW/K")
_row2("Lowest usable TES (full power, no sun)", [f"{CYC[ct]['T_TES_cutoff']:.0f}" for ct in cycles], "K")

_head("[5] OPERATION OVER THE DESIGN MONTH")
_cols()
_row2("TES at sunrise", _each(lambda r, ct: r['day']['T_TES'][0], '.0f'), "K")
_row2("TES at start of longest darkness", _each(lambda r, ct: r['night']['T_TES'][0], '.0f'), "K")
_row2("TES lowest", _each(lambda r, ct: r['cycle']['T_TES'].min(), '.0f'), "K")
_row2("Turbine inlet, min-max", _each(lambda r, ct: _rng(r['cycle']['TIT']), 's'), "K")
_row2("Cycle efficiency, min-max", _each(lambda r, ct: _rng(r['cycle']['eta'] * 100, '.1f'), 's'), "%")
_row2("Compressor inlet, min-max", _each(lambda r, ct: _rng(r['cycle']['T_comp_in']), 's'), "K")
_row2("Receiver efficiency, mean (lit)", _each(lambda r, ct: _tot(r['cycle'], 'Q_solar') /
      max(_tot(r['cycle'], 'Q_sun_in'), 1e-9) * 100, '.1f'), "%")
_row2("Solar heat curtailed (TES full)", _each(lambda r, ct: _tot(r['cycle'], 'Q_curtail') /
      max(_tot(r['cycle'], 'Q_solar'), 1e-9) * 100, '.1f'), "%")
_row2("Radiator use, peak", _each(lambda r, ct: np.nanmax(np.where(r['cycle']['running'],
      r['cycle']['rad_use'], np.nan)) * 100, '.0f'), "%")
_row2("Hours at reduced power / cycle off", _each(lambda r, ct: f"{np.sum(r['cycle']['part_load'][:-1]) * dt:.0f} / "
      f"{np.sum(~r['cycle']['running'][:-1]) * dt:.0f}", 's'), "h")
_row2("Demand served (energy)", _each(lambda r, ct: _tot(r['cycle'], 'W_elec') /
      (W_elec_target * t_lunar / 1e3) * 100, '.1f'), "%")

_head("[6] ENERGY PER DESIGN MONTH (MWh)")
_cols()
for name, fn in [
        ("Light on receiver", lambda r, ct: _tot(r['cycle'], 'Q_sun_in')),
        ("Solar heat into gas loop (net)", lambda r, ct: _tot(r['cycle'], 'Q_solar')),
        ("  curtailed", lambda r, ct: _tot(r['cycle'], 'Q_curtail')),
        ("TES heat loss (insulation)", lambda r, ct: _tot(r['cycle'], 'Q_loss')),
        ("Heat drawn by cycle", lambda r, ct: _tot(r['cycle'], 'Q_in')),
        ("  while sunlit / in longest darkness", None),
        ("Gross electricity", lambda r, ct: _tot(r['cycle'], 'W_gross')),
        ("  parasitic loads", lambda r, ct: _tot(r['cycle'], 'W_aux')),
        ("  to habitat", lambda r, ct: _tot(r['cycle'], 'W_elec')),
        ("Waste heat", lambda r, ct: _tot(r['cycle'], 'Q_out')),
        ("  to habitat / to radiator", None)]:
    if name.startswith("  while"):
        _row2(name, _each(lambda r, ct: f"{_tot(r['day'], 'Q_in'):.1f} / {_tot(r['night'], 'Q_in'):.1f}", 's'))
    elif name.startswith("  to habitat /"):
        _row2(name, _each(lambda r, ct: f"{_tot(r['cycle'], 'Q_hab'):.1f} / {_tot(r['cycle'], 'Q_rad'):.1f}", 's'))
    else:
        _row2(name, _each(fn, '.1f'))
_row2("Parasitic: controls / heliostats", _each(lambda r, ct: f"{W_aux_fixed * t_lunar / 1e3:.1f} / "
      f"{_tot(r['cycle'], 'W_helio'):.1f}", 's'))
_row2("Parasitic: blower / pump", _each(lambda r, ct: f"{_tot(r['cycle'], 'W_blower'):.1f} / "
      f"{_tot(r['cycle'], 'W_pump'):.1f}", 's'))
_row2("Parasitic share of gross output", _each(lambda r, ct: _tot(r['cycle'], 'W_aux') /
      _tot(r['cycle'], 'W_gross') * 100, '.1f'), "%")
_row2("Light on receiver -> electricity", _each(lambda r, ct: _tot(r['cycle'], 'W_elec') /
      max(_tot(r['cycle'], 'Q_sun_in'), 1e-9) * 100, '.1f'), "%")

_head("[7] HABITAT (same for both cycles)")
_row("Heating needed from the cycle, max", f"{np.max(_cy0['Q_hab_need']):.1f}", "kW")
_row("Hours the habitat needs heating", f"{np.sum(_cy0['Q_hab_need'][:-1] > 0) * dt:.0f}", "h")
_row("Habitat radiator area, max", f"{np.max(_cy0['A_hab_rad']):.0f}", "m2")
_row("Hours heat pump on", f"{np.sum(_cy0['hp_on'][:-1]) * dt:.0f}", "h")
_cols()
_row2("Habitat heating not met", _each(lambda r, ct: _tot(r['cycle'], 'Q_hab_unmet'), '.2f'), "MWh")

_head("[8] CHECKS")
_cols()
_row2("Cycle energy balance residual", _each(lambda r, ct: r['cycle_balance'], '.1e'), "kW")
_row2("TES energy balance residual", _each(lambda r, ct: r['balance'], '.1e'), "kWh")
_row2("Month repeats (TES change over month)", _each(lambda r, ct: r['periodic_error'], '.3f'), "K")
for ct in cycles:
    r = R[ct]
    w = []
    if not r['survived']:
        w.append("TES ran out")
    if abs(r['periodic_error']) > 1.0:
        w.append("not yet repeating: raise n_cycles_periodic")
    _pl = np.sum(r['cycle']['part_load'][:-1]) * dt
    if _pl > 0 or (~r['cycle']['running'][:-1]).any():
        w.append(f"reduced power {_pl:.0f} h, {_tot(r['cycle'], 'W_unserved'):.1f} MWh not served")
    if r['cycle']['Q_hab_unmet'][:-1].sum() > 1e-6:
        w.append("habitat heating not fully met")
    if abs(r['cycle_balance']) > 1e-6 or abs(r['balance']) > 1e-3:
        w.append("energy balance residual not ~0")
    rep(f"  Status ({ct}): " + ("ALL OK" if not w else "WARNING - " + "; ".join(w)))
_line()

# ##########################################################################################
# PART B -- FIGURES
# ##########################################################################################
# =====================================================================
# Plots in AIAA journal format, simple vs recuperated (run after Steps 1-11)
# AIAA "Guidelines for Journal Figures and Tables":
#   - single-column width 3.25 in; lettering >= 8 pt at that width
#   - line weight >= 0.28 pt, distinct lines (here: colour AND line style, so the
#     figures also read in grayscale), no screens lighter than 20 %
#   - no title inside the figure: the caption goes centred BELOW the figure,
#     at most 20-25 words; axis labels in words with SI units ("Temperature, K")
#   - line art at >= 600 dpi
# Each figure is drawn at 3.25 x 2.6 in, Times-style serif font. Part D puts one
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

# Time axis: the final month runs from sunrise (end of the longest darkness, Step 8).
_t_cy = R[cycles[0]]['cycle']['t']
days = (_t_cy - _t_cy[0]) / 24.0              # time since sunrise, Earth days
sunset_d = t_day / 24.0                       # start of the longest darkness
XLAB_T = 'Time since sunrise, Earth days'
# Dark steps (sun hidden at mid-step), merged into spans for shading
_dk = np.concatenate(([False], ~R[cycles[0]]['cycle']['is_day'][:-1], [False])).astype(int)
_dk0, _dk1 = np.where(np.diff(_dk) == 1)[0], np.where(np.diff(_dk) == -1)[0]
DARK_SPANS = [(days[a], days[b]) for a, b in zip(_dk0, _dk1)]
SHADE = '0.8'                                 # dark-period shading (20 % screen, AIAA minimum)

figs = []

def new_fig():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    return fig, ax

def time_axis(ax, day_only=False):
    """Time since sunrise in Earth days; dark periods shaded gray. day_only: sunrise to the
    start of the longest darkness."""
    end = sunset_d if day_only else days[-1]
    ax.set_xlim(0, end)
    ax.xaxis.set_major_locator(MultipleLocator(2 if day_only else 5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(XLAB_T)
    for a, b in DARK_SPANS:
        if a < end:
            ax.axvspan(a, min(b, end), color=SHADE, lw=0, zorder=0)

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

# 0. Sun over the design month ---------------------------------------------------------
_tt = np.linspace(_t_cy[0], _t_cy[-1], 3000)
_el, _az, _vis = sun_state(_tt)
fig, ax = new_fig()
ax.plot((_tt - _t_cy[0]) / 24, _el, color='k', ls='-', label='Sun elevation')
ax.plot((_tt - _t_cy[0]) / 24, horizon_elevation(_az), color=COL['recuperated'], ls='--',
        label='Local horizon (at sun azimuth)')
time_axis(ax)
finish(fig, ax, 'sun_elevation', 'Elevation, deg',
       f'Sun elevation and local horizon over the design month, latitude {latitude_deg:.2f} deg '
       f'(gray: sun hidden).', legend_loc='upper right', headroom=0.45)

# 1. TES temperature -------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'T_TES')
ax.axhline(T_TES_max, color=REF, lw=0.6, ls=':')
for ct in cycles:
    ax.axhline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
ax.set_ylim(min(CYC[ct]['T_TES_full_power'] for ct in cycles) - 40, T_TES_max + 20)
time_axis(ax)
finish(fig, ax, 'TES_temperature', 'TES temperature, K',
       'TES temperature over the design month (gray: dark; dotted: maximum and lowest TES giving full power).',
       legend_loc='upper left', headroom=0.15)

# 2. Stored energy -----------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'E_TES', 1e-3)
time_axis(ax)
finish(fig, ax, 'TES_energy', 'Stored energy, MWh',
       f'Usable energy stored in the TES above {T_TES_min:.0f} K over the design month (gray: dark).',
       legend_loc='upper left', headroom=0.25)

# 3. Cycle efficiency ---------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'eta', 100)
time_axis(ax)
finish(fig, ax, 'cycle_efficiency', 'Cycle efficiency, %',
       'Cycle thermal efficiency (shaft power per unit heat from the TES) over the design month.',
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
       f'Net solar heat into the gas loop and heat drawn by the cycle ({field_layout} heliostat field).',
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
           f'Parasitic electric load of the {ct}-cycle plant by source; blower and pump run whenever the cycle runs.',
           legend_loc='upper right', ncol=2)

# 6b. Gas flow in the TES loop ---------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'm_gas')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'gas_flow', 'TES gas flow, kg/s',
       f'sCO2 flow in the TES gas loop: set by the heat stored while sunlit (rise up to {dT_charge:.0f} K), by the PHX when dark.',
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
       'TES gas loop: gas leaving the PHX (back to receiver or TES) and leaving the receiver (sunlit).',
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
finish(fig, ax, 'waste_heat', 'Heat rejected, kW', 'Waste heat rejected by the plant radiator over the design month.',
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
ax.plot(days, T_surface(_cy['t'] + 0.5 * dt), color='k', lw=0.8, ls='-', label='Lunar surface')
ax.plot(days, _cy['T_sink'], color='0.45', lw=0.8, ls='-.', label='Radiator sink')
for ct in cycles:
    ax.plot(days, R[ct]['cycle']['T_comp_in'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, compressor inlet')
ax.set_ylim(0, max(float(np.max(T_surface(_cy['t'] + 0.5 * dt))), T_comp_in_max) * 1.05)
time_axis(ax)
finish(fig, ax, 'radiator_environment', 'Temperature, K',
       'Lunar surface, radiator sink (' + ('shielded' if shield else f'view factor {F_rad_surf}')
       + ') and compressor inlet temperatures.', legend_loc='upper right', headroom=0.45)

# 11. TES and turbine inlet -----------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['TIT'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, turbine inlet')
    ax.plot(days, cy['T_TES'], color=COL[ct], ls='none', marker='o', ms=1.6, markevery=max(1, int(round(24 / dt))),
            label=f'{LAB[ct]}, TES')
time_axis(ax)
finish(fig, ax, 'TES_and_TIT', 'Temperature, K',
       'Turbine inlet temperature (lines) and TES temperature (markers, every 24 h or every step).', legend_loc='lower center')

# 12. PHX gap -------------------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'dT_PHX_t')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'PHX_gap', 'TES minus turbine inlet, K',
       f'TES minus turbine inlet: {dT_TES_HX + dT_PHX:.0f} K with no sun, {dT_PHX - dT_TES_HX:.0f} K while the sun charges the TES.',
       legend_loc='upper right', headroom=0.3)

# 13. Receiver efficiency (day) -------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, np.where(cy['Q_sun_in'] > 0, cy['eta_rec_t'], np.nan) * 100,
            color=COL[ct], ls=LS[ct], label=LAB[ct])
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_efficiency', 'Receiver efficiency, %',
       f'Receiver efficiency while sunlit ({C_ratio:.0f}x concentration, emittance {eps_rec}); '
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
       'Receiver loss (reflection and thermal radiation) while sunlit.', legend_loc='best', headroom=0.15)

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
finish(fig, ax, 'habitat_radiator', 'Radiator area, m$^2$', 'Habitat radiator panel area needed over the design month.')

# 17-18. Performance versus TES temperature (no sun) --------------------------------------
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
       'Cycle efficiency versus TES temperature with no sun, radiator not limiting (dotted: full-power limit).',
       legend_loc='upper right')

fig, ax = new_fig()
for ct in cycles:
    ax.plot(*_sw(ct, 'Q_in'), color=COL[ct], ls=LS[ct], label=LAB[ct])
    ax.axvline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
ax.set_xlabel('TES temperature, K')
ax.invert_xaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
finish(fig, ax, 'Qin_vs_TES', 'Heat from TES, kW',
       f'Heat drawn from the TES for {W_elec_target:.0f} kWe with no sun, radiator not limiting (dotted: full-power limit).',
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
finish(fig, ax, 'energy_totals', 'Energy per month, MWh',
       'Energy per design month: light on receiver, solar heat collected, TES loss, heat to cycle, parasitic, habitat, radiator.',
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
say(f"Figures: {len(figs)} AIAA-format figures ({FIG_W} x {FIG_H} in, 8 pt serif)"
    + (f", saved to {FIG_DIR}/ as PDF and 600-dpi PNG" if SAVE_FIGS else ''))

# ##########################################################################################
# PART C -- REFERENCES
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
    ("R20", "latitude_deg, horizon (Step 2), illumination check", "89.45 S; 88.1 % lit, 4.6 d longest dark",
     "Glaser, P., Scholten, F., De Rosa, D., Marco Figuera, R., Oberst, J., Mazarico, E., Neumann, G.A., "
     "Robinson, M.S. (2014). Illumination conditions at the lunar south pole using high resolution "
     "Digital Terrain Models from LOLA. Icarus 243, 78-90. doi:10.1016/j.icarus.2014.08.013"),
    ("R21", "horizon method (Step 8)", "horizon elevation vs azimuth",
     "Mazarico, E., Neumann, G.A., Smith, D.E., Zuber, M.T., Torrence, M.H. (2011). Illumination "
     "conditions of the lunar polar regions using LOLA topography. Icarus 211, 1066-1081. "
     "doi:10.1016/j.icarus.2010.10.030"),
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
    ("Synthetic horizon (dip 2.0 deg, one 4.0 deg range)", "calibrated to Connecting Ridge illumination [R20]; a LOLA horizon can replace it"),
    ("Design month = darkest synodic month, repeated", "sizing case: the longest darkness of the year every month (conservative)"),
    ("heliostat_spacing_ratio = None", "no mutual shading of heliostats at low sun (optimistic at the pole)"),
    ("T_ground = 250 K (equatorial)", "polar subsurface is colder: TES heat loss is underestimated"),
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

say(f"References: {len(REFERENCES)} references, {len(ASSUMED)} assumptions and "
    f"{len(VALIDATION)} validation sources listed for the PDF")

# ##########################################################################################
# PART D -- RESULTS PDF
# ##########################################################################################
# =====================================================================
# ONE RESULTS PDF: summary page, report, figures, references (run log only if PDF_RUN_LOG)
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
    s = [f"Run {datetime.datetime.now():%Y-%m-%d %H:%M} | site {SITE.replace('_', ' ')}, latitude "
         f"{latitude_deg:.2f} deg | time step {dt:g} h",
         f"Sun: lit {100*SITE_STATS['illum_fraction']:.1f} % of the year | design month {t_lunar:.0f} h: "
         f"sunlit {t_day:.0f} h, longest darkness {t_night:.0f} h",
         f"Demand {W_elec_target:.0f} kWe | TES {T_TES_max:.0f} -> {T_TES_min:.0f} K | "
         f"compressor inlet {T_comp_in:.0f}-{T_comp_in_max:.0f} K | {P_low/1e6:.0f} -> {P_high/1e6:.0f} MPa",
         "", f"{'':40}" + "".join(f"{ct.capitalize():>16}" for ct in cycles)]
    def row(name, fn, fmt, unit=''):
        s.append(f"{name:<40}" + "".join(f"{format(fn(ct), fmt):>16}" for ct in cycles) + f"  {unit}")
    row("Cycle efficiency at full charge", lambda ct: cycle_at(TIT_max, ct)['eta'] * 100, '.1f', '%')
    row("TES mass", lambda ct: R[ct]['m_TES'] / 1e3, ',.0f', 't')
    row("TES volume", lambda ct: R[ct]['V_TES'], ',.0f', 'm3')
    row("Collector (mirror) area", lambda ct: R[ct]['A_col'], ',.0f', 'm2')
    row("Radiator panel area", lambda ct: A_RAD[ct], ',.0f', 'm2')
    row("Radiator mass", lambda ct: A_RAD[ct] * rad_areal_mass / 1e3, '.1f', 't')
    row("TES lowest in operation", lambda ct: R[ct]['cycle']['T_TES'].min(), '.0f', 'K')
    row("Parasitic share of gross output", lambda ct: np.sum(R[ct]['cycle']['W_aux'][:-1]) /
        np.sum(R[ct]['cycle']['W_gross'][:-1]) * 100, '.1f', '%')
    row("Demand served (energy)", lambda ct: np.sum(R[ct]['cycle']['W_elec'][:-1]) * dt /
        (W_elec_target * t_lunar) * 100, '.1f', '%')
    row("Light on receiver -> electricity", lambda ct: np.sum(R[ct]['cycle']['W_elec'][:-1]) /
        np.sum(R[ct]['cycle']['Q_sun_in'][:-1]) * 100, '.1f', '%')
    s += [""] + [l.strip() for l in REPORT if l.strip().startswith('Status')]
    return s

pages = 0
with PdfPages(OUT_PDF) as pdf:
    pages += _text_pages(pdf, 'Lunar solar-TES sCO2 power system: summary', _summary_lines(),
                         fontsize=8.0, width=110)
    pages += _text_pages(pdf, 'Report', REPORT)
    if PDF_RUN_LOG:
        pages += _text_pages(pdf, 'Run log (messages from Steps 2-11)', RUN_LOG)
    for k, (name, f, caption) in enumerate(figs, 1):
        # one figure per page: the AIAA-size figure (3.25 in) shown at 2x, caption centred below
        img = plt.imread(io.BytesIO(FIG_PNG[k - 1]))          # rendered once in Part B
        page = plt.figure(figsize=(8.5, 11))
        w_in = 2 * FIG_W
        h_in = w_in * img.shape[0] / img.shape[1]
        ax_img = page.add_axes([(8.5 - w_in) / 2 / 8.5, 1 - (1.0 + h_in) / 11, w_in / 8.5, h_in / 11])
        ax_img.imshow(img, interpolation='none')        # stored as is: no resampling (fast, sharp)
        ax_img.axis('off')
        page.text(0.5, 1 - (1.25 + h_in) / 11, f"Fig. {k}  {caption}", ha='center', va='top',
                  fontsize=11, family='serif', weight='bold', wrap=True)
        pdf.savefig(page)
        plt.close(page)
        pages += 1
    if PDF_REFERENCES:
        pages += _text_pages(pdf, 'References, assumptions and validation', REFS_LINES)
    info = pdf.infodict()
    info['Title'] = 'Lunar solar-TES sCO2 power system results'
    info['Author'] = 'Lunar_Regolith_TES_Cycle notebook'

print(f"Saved {OUT_PDF}: {pages} pages (summary, report, {len(figs)} figures"
      + (", run log" if PDF_RUN_LOG else "") + (", references" if PDF_REFERENCES else "") + ")")