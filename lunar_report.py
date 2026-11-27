# lunar_report.py -- report, figures, references and the results PDF, in one go
# Loaded by Step 12 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses the results of Steps 1-11 (R, CYC, A_RAD, ...).
#   Part A  text report (short: the numbers for the paper)
#   Part B  AIAA figures, saved to FIG_DIR
#   Part C  references, assumptions, validation
#   Part D  one results PDF: summary page, report, figures, references
# Site: lunar south pole (Connecting Ridge). Time in the report and the figures = time since
# sunrise (the end of the longest darkness of the design month, Step 8). "Dark" = sun hidden
# at the middle of the step.
# TES: N_segments brick segments, END A hot (gas to the PHX), END B cold (gas to the receiver).

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
    a = np.asarray(a, dtype=float)
    return f"{np.nanmin(a):{fmt}}-{np.nanmax(a):{fmt}}" if np.isfinite(a).any() else '-'
def _pct(a, b):
    return 100 * a / max(b, 1e-9)

_cy0 = R[cycles[0]]['cycle']

_line()
rep("  LUNAR POWER SYSTEM REPORT")
rep("  Solar tower -> sintered-regolith TES (two gas loops) -> sCO2 Brayton (simple vs recuperated)")
rep(f"  Site: lunar south pole, Connecting Ridge, latitude {latitude_deg:.2f} deg | time step {dt:g} h")
_line()

_head("[1] SITE AND SUN")
_row("Sunlit over one draconic year", f"{100*SITE_STATS['illum_fraction']:.1f}", "%")
_row("Longest darkness over the year", f"{SITE_STATS['longest_dark_h']:.0f}", "h")
_row("Sun elevation, min / max", f"{SITE_STATS['el_min_deg']:.1f} / {SITE_STATS['el_max_deg']:.1f}", "deg")
_row("Design month (darkest of the year)", f"{t_lunar:.0f}", "h")
_row("  sunlit / longest darkness", f"{t_day:.0f} / {t_night:.0f}", "h")
_vis_m = sun_state(_cy0['t'][:-1] + 0.5 * dt)[2]
_row("  full sun / sun partly hidden / dark", f"{np.sum(_vis_m >= 0.999) * dt:.0f} / "
     f"{np.sum((_vis_m > 0) & (_vis_m < 0.999)) * dt:.0f} / {np.sum(_vis_m <= 0) * dt:.0f}", "h")
_row("Mean heliostat cosine factor (lit)", f"{SITE_STATS['mean_cosine_lit']:.3f}")
_row("Local horizon", "synthetic (2.0 deg dip, 4.0 deg range)")

_head("[2] MAIN INPUTS (all inputs: lunar_inputs.py)")
_row("Electrical demand / generator efficiency", f"{W_elec_target:.0f} kWe / {eta_gen:.2f}")
_row("TES bricks: range / regolith / density", f"{T_TES_max:.0f}-{T_TES_min:.0f} K / {regolith_type} / {rho_TES:.0f} kg/m3")
_row("TES segments / dp per segment / mass margin", f"{N_segments} / {dp_seg/1e5:g} bar / {TES_mass_margin:.2f}")
_row("TES gas loops: pressure / receiver outlet", f"{P_loop/1e6:.0f} MPa / {T_rec_set:.0f} K")
_row("Smallest receiver rise / blower efficiency", f"{dT_charge:.0f} K / {eta_blower:.2f}")
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
_row2("CO2 flow (cycle) / gas flow (blower 2)", _each(lambda r, ct: f"{r['p_design']['m_dot']:.2f} / "
      f"{r['p_design']['m_gas']:.2f}", 's'), "kg/s")
_row2("Gas back to TES end B", _each(lambda r, ct: r['p_design']['T_gas_cold'], '.0f'), "K")
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
_row2("TES volume / per segment", _each(lambda r, ct: f"{r['V_TES']:,.0f} / {r['V_TES']/N_segments:,.0f}", 's'), "m3")
_row2("Collector (mirror) area", _each(lambda r, ct: r['A_col'], ',.0f'), "m2")
_row2("Radiator panel area (peak use)", [f"{A_RAD[ct]:,.0f}" for ct in cycles], "m2")
_row2("Radiator mass", [f"{A_RAD[ct] * rad_areal_mass / 1e3:.1f}" for ct in cycles], "t")
_row2("PHX size UA", [f"{UA_PHX[ct] / 1e3:.1f}" for ct in cycles], "kW/K")
_row2("Store floor: lowest gas from end A", [f"{CYC[ct]['T_TES_cutoff']:.0f}" for ct in cycles], "K")
_row2("Gas loop pressure drop (each loop)", [f"{loop_dp()/1e5:.2f}"] * len(cycles), "bar")

_head("[5] OPERATION OVER THE DESIGN MONTH")
_cols()
_row2("Gas from end A at sunrise", _each(lambda r, ct: r['day']['T_A'][0], '.0f'), "K")
_row2("Gas from end A at start of longest darkness", _each(lambda r, ct: r['night']['T_A'][0], '.0f'), "K")
_row2("Gas from end A, lowest", _each(lambda r, ct: r['cycle']['T_A'].min(), '.0f'), "K")
_row2("Gas from end B (to receiver), min-max", _each(lambda r, ct: _rng(r['cycle']['T_B']), 's'), "K")
_row2("Store charge at start of longest darkness", _each(lambda r, ct: _pct(r['night']['E_TES'][0],
      stored_energy(np.full(N_segments, T_TES_max), tes_geometry(r['m_TES']))), '.0f'), "%")
_row2("Turbine inlet, min-max", _each(lambda r, ct: _rng(r['cycle']['TIT']), 's'), "K")
_row2("Cycle efficiency, min-max", _each(lambda r, ct: _rng(r['cycle']['eta'] * 100, '.1f'), 's'), "%")
_row2("Compressor inlet, min-max", _each(lambda r, ct: _rng(r['cycle']['T_comp_in']), 's'), "K")
_row2("Receiver efficiency, mean (lit)", _each(lambda r, ct: _pct(_tot(r['cycle'], 'Q_solar'),
      _tot(r['cycle'], 'Q_sun_in')), '.1f'), "%")
_row2("Solar heat curtailed (store full)", _each(lambda r, ct: _pct(_tot(r['cycle'], 'Q_curtail'),
      _tot(r['cycle'], 'Q_solar')), '.1f'), "%")
_row2("Blower 1 (charging) peak", _each(lambda r, ct: np.nanmax(r['cycle']['W_bl1']), '.1f'), "kWe")
_row2("Blower 2 (discharging) peak", _each(lambda r, ct: np.nanmax(r['cycle']['W_blower']), '.1f'), "kWe")
_row2("Radiator use, peak", _each(lambda r, ct: np.nanmax(np.where(r['cycle']['running'],
      r['cycle']['rad_use'], np.nan)) * 100, '.0f'), "%")
_row2("Hours at reduced power / cycle off", _each(lambda r, ct: f"{np.sum(r['cycle']['part_load'][:-1]) * dt:.0f} / "
      f"{np.sum(~r['cycle']['running'][:-1]) * dt:.0f}", 's'), "h")
_row2("Demand served (energy)", _each(lambda r, ct: _pct(_tot(r['cycle'], 'W_elec'),
      W_elec_target * t_lunar / 1e3), '.1f'), "%")

_head("[6] ENERGY PER DESIGN MONTH (MWh)")
_cols()
_row2("Light on receiver", _each(lambda r, ct: _tot(r['cycle'], 'Q_sun_in'), '.1f'))
_row2("Solar heat into the charging gas (net)", _each(lambda r, ct: _tot(r['cycle'], 'Q_solar'), '.1f'))
_row2("  stored in the bricks / curtailed", _each(lambda r, ct: f"{_tot(r['cycle'], 'Q_ch'):.1f} / "
      f"{_tot(r['cycle'], 'Q_curtail'):.1f}", 's'))
_row2("Taken from the bricks (discharging gas)", _each(lambda r, ct: _tot(r['cycle'], 'Q_dis'), '.1f'))
_row2("TES heat loss (insulation)", _each(lambda r, ct: _tot(r['cycle'], 'Q_loss'), '.1f'))
_row2("Heat into the cycle (PHX)", _each(lambda r, ct: _tot(r['cycle'], 'Q_in'), '.1f'))
_row2("  while sunlit / in longest darkness", _each(lambda r, ct: f"{_tot(r['day'], 'Q_in'):.1f} / "
      f"{_tot(r['night'], 'Q_in'):.1f}", 's'))
_row2("Gross electricity", _each(lambda r, ct: _tot(r['cycle'], 'W_gross'), '.1f'))
_row2("  parasitic loads / to habitat", _each(lambda r, ct: f"{_tot(r['cycle'], 'W_aux'):.1f} / "
      f"{_tot(r['cycle'], 'W_elec'):.1f}", 's'))
_row2("Waste heat", _each(lambda r, ct: _tot(r['cycle'], 'Q_out'), '.1f'))
_row2("  to habitat / to radiator", _each(lambda r, ct: f"{_tot(r['cycle'], 'Q_hab'):.1f} / "
      f"{_tot(r['cycle'], 'Q_rad'):.1f}", 's'))
_row2("Parasitic: controls / heliostats", _each(lambda r, ct: f"{W_aux_fixed * t_lunar / 1e3:.1f} / "
      f"{_tot(r['cycle'], 'W_helio'):.1f}", 's'))
_row2("Parasitic: blower 1 / blower 2", _each(lambda r, ct: f"{_tot(r['cycle'], 'W_bl1'):.1f} / "
      f"{_tot(r['cycle'], 'W_blower'):.1f}", 's'))
_row2("Parasitic: coolant pump", _each(lambda r, ct: _tot(r['cycle'], 'W_pump'), '.1f'))
_row2("Parasitic share of gross output", _each(lambda r, ct: _pct(_tot(r['cycle'], 'W_aux'),
      _tot(r['cycle'], 'W_gross')), '.1f'), "%")
_row2("Light on receiver -> electricity", _each(lambda r, ct: _pct(_tot(r['cycle'], 'W_elec'),
      _tot(r['cycle'], 'Q_sun_in')), '.1f'), "%")

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
_row2("Sun balance residual", _each(lambda r, ct: r['sun_balance'], '.1e'), "kWh")
_row2("Heat from bricks vs heat to cycle", _each(lambda r, ct: r['phx_diff'], '+.2f'), "%")
_row2("Month repeats (stored energy change)", _each(lambda r, ct: r['periodic_error'], '.1f'), "kWh")
for ct in cycles:
    r = R[ct]
    w = []
    if not r['survived']:
        w.append("TES ran out")
    if abs(r['periodic_error']) > 1e-3 * max(np.max(r['cycle']['E_TES']), 1.0):
        w.append("not yet repeating: raise n_cycles_periodic")
    _pl = np.sum(r['cycle']['part_load'][:-1]) * dt
    if _pl > 0 or (~r['cycle']['running'][:-1]).any():
        w.append(f"reduced power {_pl:.0f} h, {_tot(r['cycle'], 'W_unserved'):.1f} MWh not served")
    if r['cycle']['Q_hab_unmet'][:-1].sum() > 1e-6:
        w.append("habitat heating not fully met")
    if abs(r['cycle_balance']) > 1e-6 or abs(r['balance']) > 1e-5 * max(_tot(r['cycle'], 'Q_ch') * 1e3, 1.0):
        w.append("energy balance residual not ~0")
    rep(f"  Status ({ct}): " + ("ALL OK" if not w else "WARNING - " + "; ".join(w)))
_line()

# ##########################################################################################
# PART B -- FIGURES
# ##########################################################################################
# AIAA "Guidelines for Journal Figures and Tables":
#   - single-column width 3.25 in; lettering >= 8 pt at that width
#   - line weight >= 0.28 pt, distinct lines (colour AND line style, so the figures also
#     read in grayscale), no screens lighter than 20 %
#   - no title inside the figure: the caption goes centred BELOW the figure, short;
#     axis labels in words with SI units ("Temperature, K")
#   - line art at >= 600 dpi
# Each figure is drawn at 3.25 x 2.6 in, Times-style serif font. Part D puts one figure per
# PDF page with its caption, and (SAVE_FIGS) also writes every figure at true size to FIG_DIR
# as vector PDF and 600-dpi PNG for the paper.  figs = list of (name, figure, caption)
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
import os

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
LS2 = {'simple': '-.', 'recuperated': ':'}    # second line of the same cycle
LAB = {'simple': 'Simple', 'recuperated': 'Recuperated'}
REF = '0.35'                                  # reference lines: dark gray (>= 20 % screen)

# Time axis: the final month runs from sunrise (end of the longest darkness, Step 8).
_t_cy = _cy0['t']
days = (_t_cy - _t_cy[0]) / 24.0              # time since sunrise, Earth days
sunset_d = t_day / 24.0                       # start of the longest darkness
XLAB_T = 'Time since sunrise, Earth days'
# Sun at the middle of each step: 0 = full sun (white), 1 = sun partly hidden behind the
# horizon (light gray), 2 = dark (medium gray). AIAA: no screen lighter than 20 %.
_vis_step = sun_state(_cy0['t'][:-1] + 0.5 * dt)[2]
SUN_CLASS = np.where(_vis_step <= 0, 2, np.where(_vis_step < 0.999, 1, 0))
def _spans(mask):
    m = np.concatenate(([0], mask.astype(int), [0]))
    a, b = np.where(np.diff(m) == 1)[0], np.where(np.diff(m) == -1)[0]
    return [(days[i], days[j]) for i, j in zip(a, b)]
SHADES = [(_spans(SUN_CLASS == 1), '0.8'),    # partly sunlit: light gray (20 % screen)
          (_spans(SUN_CLASS == 2), '0.6')]    # dark: medium gray (40 % screen)
SHADE_NOTE = ' Light gray: sun partly hidden; gray: dark.'

figs = []

def new_fig():
    return plt.subplots(figsize=(FIG_W, FIG_H))

def time_axis(ax, day_only=False):
    """Time since sunrise in Earth days; dark periods shaded gray. day_only: sunrise to the
    start of the longest darkness."""
    end = sunset_d if day_only else days[-1]
    ax.set_xlim(0, end)
    ax.xaxis.set_major_locator(MultipleLocator(2 if day_only else 5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(XLAB_T)
    for spans, shade in SHADES:
        for a, b in spans:
            if a < end:
                ax.axvspan(a, min(b, end), color=shade, lw=0, zorder=0)
    ax._sun_shaded = True

def both(ax, key, scale=1.0, day_mask=False, **kw):
    for ct in cycles:
        cy = R[ct]['cycle']
        y = cy[key] * scale
        if day_mask:
            y = np.where(cy['Q_sun_in'] > 0, y, np.nan)
        ax.plot(days, y, color=COL[ct], ls=LS[ct], label=LAB[ct], **kw)

def finish(fig, ax, name, ylabel, caption, legend_loc='best', ncol=1, headroom=0.0):
    ax.set_ylabel(ylabel)
    if headroom:                       # space above the data so the legend sits clear of the lines
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + (hi - lo) * headroom)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc=legend_loc, ncol=ncol)
    fig.tight_layout(pad=0.3)
    if getattr(ax, '_sun_shaded', False):
        caption += SHADE_NOTE
    figs.append((name, fig, caption))

# --- Sun ----------------------------------------------------------------------------------
_tt = np.linspace(_t_cy[0], _t_cy[-1], 3000)
_el, _az, _vis = sun_state(_tt)
fig, ax = new_fig()
ax.plot((_tt - _t_cy[0]) / 24, _el, color='k', ls='-', label='Sun elevation')
ax.plot((_tt - _t_cy[0]) / 24, horizon_elevation(_az), color=COL['recuperated'], ls='--',
        label='Local horizon (at sun azimuth)')
time_axis(ax)
finish(fig, ax, 'sun_elevation', 'Elevation, deg',
       f'Sun elevation and local horizon over the design month, latitude {latitude_deg:.2f} deg.', legend_loc='upper right', headroom=0.45)

# --- Share of the sun disk above the horizon -------------------------------------------------
fig, ax = new_fig()
ax.plot((_tt - _t_cy[0]) / 24, _vis * 100, color='k', ls='-')
ax.axhline(100, color=REF, lw=0.6, ls=':')
ax.set_ylim(0, 105)
time_axis(ax)
finish(fig, ax, 'sun_visible', 'Sun disk visible, %',
       'Share of the sun disk above the local horizon over the design month.')

# --- TES: gas from end A and end B ------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['T_A'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, end A')
    ax.plot(days, cy['T_B'], color=COL[ct], ls=LS2[ct], label=f'{LAB[ct]}, end B')
ax.axhline(T_TES_max, color=REF, lw=0.6, ls=':')
time_axis(ax)
finish(fig, ax, 'TES_end_temperatures', 'Gas temperature, K',
       'Gas leaving the two ends of the TES: end A (hot, to the PHX) and end B (cold, to the receiver).',
       legend_loc='upper center', ncol=2, headroom=0.5)

# --- TES: temperature along the segments, one figure per cycle -------------------------------
_seg = np.arange(1, N_segments + 1)
for ct in cycles:
    r = R[ct]
    fig, ax = new_fig()
    _nt = r['night']['T_seg']                          # its last point = next sunrise
    for lbl, T, ls in [('Start of longest darkness', _nt[0], '--'),
                       ('Middle of longest darkness', _nt[len(_nt) // 2], ':'),
                       ('Sunrise (end of darkness)', _nt[-1], '-')]:
        ax.plot(_seg, T, color=COL[ct], ls=ls, marker='o', ms=2.5, label=lbl)
    ax.axhline(CYC[ct]['T_TES_cutoff'], color=REF, lw=0.6, ls='-.')
    ax.set_xlabel('Segment (1 = end A, hot)')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_minor_locator(plt.NullLocator())
    finish(fig, ax, f'TES_profile_{ct}', 'Brick temperature, K',
           f'Brick temperature along the TES, {ct} cycle (dash-dot: store floor for end A).',
           legend_loc='upper right', headroom=0.4)

# --- Stored energy ----------------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'E_TES', 1e-3)
time_axis(ax)
finish(fig, ax, 'TES_energy', 'Stored energy, MWh',
       f'Energy stored in the bricks above {T_TES_min:.0f} K over the design month.',
       legend_loc='upper left', headroom=0.25)

# --- Cycle efficiency ---------------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'eta', 100)
time_axis(ax)
finish(fig, ax, 'cycle_efficiency', 'Cycle efficiency, %',
       'Cycle thermal efficiency (shaft power per unit heat into the cycle) over the design month.',
       legend_loc='upper left', headroom=0.3)

# --- Turbine inlet and gas from end A ----------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['TIT'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, turbine inlet')
    ax.plot(days, cy['T_A'], color=COL[ct], ls='none', marker='o', ms=1.6,
            markevery=max(1, int(round(24 / dt))), label=f'{LAB[ct]}, gas from end A')
time_axis(ax)
finish(fig, ax, 'TES_and_TIT', 'Temperature, K',
       f'Turbine inlet (lines) and gas from TES end A (markers): {dT_PHX:.0f} K PHX approach.',
       legend_loc='lower center')

# --- Heat flows, one figure per cycle ----------------------------------------------------------
for ct in cycles:
    cy = R[ct]['cycle']
    fig, ax = new_fig()
    ax.plot(days, cy['Q_solar'], color='k', ls='-', label='Solar heat (net)')
    ax.plot(days, cy['Q_ch'], color=COL[ct], ls='--', label='Stored in bricks')
    ax.plot(days, cy['Q_curtail'], color='0.45', ls=':', label='Curtailed')
    ax.plot(days, cy['Q_in'], color=COL[ct], ls='-.', label='Into the cycle')
    ax.set_ylim(bottom=0)
    time_axis(ax)
    finish(fig, ax, f'heat_flows_{ct}', 'Thermal power, kW',
           f'Heat flows of the {ct}-cycle plant: solar heat, stored in the bricks, curtailed, into the cycle.',
           legend_loc='upper right', ncol=2, headroom=0.35)

# --- Parasitic loads, one figure per cycle -------------------------------------------------------
_ymax = max(np.nanmax(R[ct]['cycle']['W_aux']) for ct in cycles) * 1.45
for ct in cycles:
    cy = R[ct]['cycle']
    fig, ax = new_fig()
    ax.plot(days, cy['W_aux'], color='k', lw=1.2, ls='-', label='Total')
    ax.plot(days, cy['W_bl1'], color=COL[ct], ls='--', label='Blower 1 (charging)')
    ax.plot(days, cy['W_blower'], color=COL[ct], ls='-.', label='Blower 2 (discharging)')
    ax.plot(days, cy['W_helio'], color='0.35', ls=(0, (6, 2, 1, 2, 1, 2)), lw=0.9, label='Heliostat drives')
    ax.plot(days, np.full_like(days, W_aux_fixed), color='0.35', ls=':', label='Controls')
    ax.plot(days, cy['W_pump'], color='0.55', ls='-', lw=0.8, label='Coolant pump')
    ax.set_ylim(0, _ymax)
    time_axis(ax)
    finish(fig, ax, f'parasitic_{ct}', 'Parasitic load, kWe',
           f'Parasitic electric load of the {ct}-cycle plant by source.',
           legend_loc='upper right', ncol=2)

# --- Gas flows in the two TES loops -------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['m_ch'], color=COL[ct], ls=LS2[ct], label=f'{LAB[ct]}, blower 1')
    ax.plot(days, np.where(cy['running'], cy['m_gas'], 0.0), color=COL[ct], ls=LS[ct],
            label=f'{LAB[ct]}, blower 2')
ax.set_ylim(bottom=0)
time_axis(ax)
finish(fig, ax, 'gas_flows', 'Gas flow, kg/s',
       f'sCO2 flow in the charging loop (blower 1, receiver outlet held at {T_rec_set:.0f} K) '
       f'and the discharging loop (blower 2, set by the cycle).', legend_loc='upper right', ncol=2,
       headroom=0.6)

# --- Coolant-loop temperatures --------------------------------------------------------------------
fig, ax = new_fig()
for ct in cycles:
    cy = R[ct]['cycle']
    ax.plot(days, cy['T_w_hot'], color=COL[ct], ls=LS[ct], label=f'{LAB[ct]}, from precooler')
    ax.plot(days, cy['T_w_cold'], color=COL[ct], ls=LS2[ct], label=f'{LAB[ct]}, to precooler')
ax.axhline(T_cool_max, color=REF, lw=0.6, ls=':')
time_axis(ax)
finish(fig, ax, 'coolant_temperature', 'Water temperature, K',
       f'Coolant-loop water leaving and entering the precooler (dotted: {T_cool_max:.0f} K limit).',
       legend_loc='center right', headroom=0.1)

# --- Compressor inlet -----------------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'T_comp_in')
ax.axhline(T_comp_in_max, color=REF, lw=0.6, ls=':')
ax.set_ylim(T_comp_in - 5, T_comp_in_max + 5)
time_axis(ax)
finish(fig, ax, 'compressor_inlet', 'Compressor inlet temperature, K',
       f'Compressor inlet temperature set by the radiator (allowed {T_comp_in:.0f}-{T_comp_in_max:.0f} K, dotted: upper limit).',
       legend_loc='center right')

# --- Waste heat to radiator -----------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'Q_rad')
time_axis(ax)
finish(fig, ax, 'waste_heat', 'Heat rejected, kW', 'Waste heat rejected by the plant radiator over the design month.',
       legend_loc='upper center', ncol=2, headroom=0.25)

# --- Radiator use -------------------------------------------------------------------------------
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

# --- Receiver (sunlit only) -----------------------------------------------------------------------
fig, ax = new_fig()
both(ax, 'eta_rec_t', 100, day_mask=True)
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_efficiency', 'Receiver efficiency, %',
       f'Receiver efficiency while sunlit ({C_ratio:.0f}x concentration, emittance {eps_rec}).',
       legend_loc='best')

fig, ax = new_fig()
both(ax, 'T_rec', day_mask=True)
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_temperature', 'Receiver temperature, K',
       f'Receiver wall temperature: mean of gas in (end B) and out ({T_rec_set:.0f} K) plus {dT_wall:.0f} K.',
       legend_loc='best')

fig, ax = new_fig()
both(ax, 'Q_rec_loss', day_mask=True)
ax.set_ylim(bottom=0)
time_axis(ax, day_only=True)
finish(fig, ax, 'receiver_loss', 'Receiver loss, kW',
       'Receiver loss (reflection and thermal radiation) while sunlit.', legend_loc='best', headroom=0.15)

# --- Habitat --------------------------------------------------------------------------------------
fig, ax = new_fig()
ax.plot(days, _cy0['Q_hab_cool'], color='k', ls='-', label='Cooling load')
ax.plot(days, _cy0['Q_hab_wall'], color='0.45', ls='--', label='Wall loss (+) / gain (-)')
ax.plot(days, _cy0['Q_hab_need'], color=COL['recuperated'], ls=':', label='Heating from cycle')
ax.plot(days, _cy0['W_hp'], color=COL['simple'], ls='-.', label='Heat-pump power')
ax.axhline(0, color=REF, lw=0.4)
time_axis(ax)
finish(fig, ax, 'habitat_balance', 'Power, kW',
       f'Habitat heat balance ({n_crew} crew, {V_hab:.0f} m$^3$, {T_hab:.0f} K); the same for both cycles.',
       legend_loc='center right')

# --- Performance versus gas from end A (no sun) -----------------------------------------------------
# Swept with the radiator NOT limiting (compressor inlet held at T_comp_in), so the curves show
# the plant itself over the whole store; dotted lines mark the lowest end-A gas at which the
# installed radiator still gives full power.
_sweep, _A_keep = {}, dict(A_RAD)
for ct in cycles:
    A_RAD[ct] = np.inf
    T_sw = np.linspace(CYC[ct]['T_TES_cutoff'], T_TES_max, 120)
    P_sw = [plant_point(T, ct=ct) for T in T_sw]
    _sweep[ct] = (T_sw, [p if (p['running'] and not p['part_load']) else None for p in P_sw])
A_RAD.update(_A_keep)
def _sw(ct, key, scale=1.0):
    T_sw, P_sw = _sweep[ct]
    return T_sw, np.array([p[key] * scale if p is not None else np.nan for p in P_sw])
for key, scale, name, ylab, cap, loc in [
        ('eta', 100, 'eta_vs_TES', 'Cycle efficiency, %', 'Cycle efficiency', 'upper right'),
        ('Q_in', 1, 'Qin_vs_TES', 'Heat into the cycle, kW',
         f'Heat into the cycle for {W_elec_target:.0f} kWe', 'upper left'),
        ('W_blower', 1, 'blower2_vs_TES', 'Blower 2 power, kWe', 'Discharging blower power', 'upper left'),
        ('T_gas_cold', 1, 'TB_return_vs_TES', 'Gas back to end B, K', 'Gas returned to TES end B', 'upper right')]:
    fig, ax = new_fig()
    for ct in cycles:
        ax.plot(*_sw(ct, key, scale), color=COL[ct], ls=LS[ct], label=LAB[ct])
        ax.axvline(CYC[ct]['T_TES_full_power'], color=COL[ct], lw=0.6, ls=':')
    ax.set_xlabel('Gas from TES end A, K')
    ax.invert_xaxis()
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    finish(fig, ax, name, ylab,
           f'{cap} versus gas from TES end A, no sun, radiator not limiting (dotted: full-power limit).',
           legend_loc=loc)

# --- T-s diagram ----------------------------------------------------------------------------------
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

# --- Energy totals --------------------------------------------------------------------------------
items = [('Q_sun_in', 'Light'), ('Q_solar', 'Solar\nheat'), ('Q_ch', 'Stored'),
         ('Q_loss', 'TES\nloss'), ('Q_in', 'To\ncycle'), ('W_aux', 'Para-\nsitic'),
         ('W_elec', 'Habi-\ntat'), ('Q_rad', 'Radi-\nator')]
x = np.arange(len(items)); w = 0.38
fig, ax = new_fig()
_hatch = {'simple': '', 'recuperated': '////'}
for i, ct in enumerate(cycles):
    off = (-1 if i == 0 else 1) * (w / 2 + 0.01)
    ax.bar(x + off, [_tot(R[ct]['cycle'], k) for k, _ in items], w, color='white' if i else COL[ct],
           edgecolor=COL[ct], hatch=_hatch[ct], lw=0.6, label=LAB[ct])
ax.set_xticks(x, [lbl for _, lbl in items])
ax.tick_params(axis='x', which='both', top=False, bottom=False)
ax.xaxis.set_minor_locator(plt.NullLocator())
ax.set_ylim(bottom=0)
finish(fig, ax, 'energy_totals', 'Energy per month, MWh',
       'Energy per design month, from the light on the receiver to the habitat and the radiator.',
       legend_loc='upper right')

# --- Render each figure ONCE to a 600-dpi PNG; the same image goes to FIG_DIR and the PDF ---
FIG_PNG = []
if SAVE_FIGS:
    os.makedirs(FIG_DIR, exist_ok=True)
    for _old in os.listdir(FIG_DIR):                  # remove figures of earlier runs (numbers shift)
        if _old.startswith('Fig') and _old.endswith(('.pdf', '.png')):
            os.remove(os.path.join(FIG_DIR, _old))
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
# Each input is tagged in lunar_inputs.py with its reference number [R..] or [A] (assumed:
# no single source; state it in the paper and cover it with a sensitivity study).
REFERENCES = [
    # tag,  inputs,  value used,  source
    ("R1", "G_sun", "1361 W/m2",
     "Kopp, G., Lean, J.L. (2011). A new, lower value of total solar irradiance: Evidence and "
     "climate significance. Geophysical Research Letters 38, L01706. doi:10.1029/2010GL045777"),
    ("R2", "rho_mirror", "0.95",
     "Typical silvered-glass heliostat reflectance (~95 %); see e.g. Stine, W.B., Geyer, M. "
     "(2001). Power From The Sun, Ch. 10 (central receivers)"),
    ("R3", "f_dust", "0.93",
     "Terrestrial heliostat soiling losses of 6.3-8.8 % (Strachan et al., 1993, Sandia). Lunar "
     "dust loss is not known: assumed"),
    ("R4", "alpha_rec, eps_rec", "0.95, 0.85",
     "Ho, C.K., Iverson, B.D. (2014). Review of high-temperature central receiver designs for "
     "concentrating solar power. Renewable and Sustainable Energy Reviews 29, 835-846. "
     "doi:10.1016/j.rser.2013.08.099"),
    ("R5", "C_ratio", "500 (towers 500-1000)",
     "Ho & Iverson (2014) [R4]; Ho, C.K. (2017). Advances in central receivers for concentrating "
     "solar applications. Solar Energy 152, 38-56. doi:10.1016/j.solener.2017.03.048"),
    ("R6", "cosine factor (Step 8)", "cos(0.5 * angle sun-receiver)",
     "Heliostat normal bisects the sun and receiver directions: Stine & Geyer (2001) [R2], Ch. 10"),
    ("R7", "Regolith cp(T), k(T)", "Schreiner correlations",
     "Paper Table 3, ref. [16]: Schreiner, S.S., Dominguez, J.A., Sibille, L., Hoffman, J.A. "
     "(2016). Thermophysical property models for lunar regolith. Advances in Space Research "
     "57(5), 1209-1222"),
    ("R8", "t_lunar, T_surf_noon", "29.53 d, 392 K",
     "Williams, J.-P. et al. (2017). The global surface temperatures of the Moon as measured by "
     "Diviner. Icarus 283, 300-325. doi:10.1016/j.icarus.2016.08.012; Vasavada, A.R. et al. "
     "(2012). JGR Planets 117, E00H18. doi:10.1029/2011JE003987"),
    ("R10", "eta_comp, eta_turb", "0.85, 0.90",
     "NETL (2014). An Assessment of Supercritical CO2 Power Cycles Integrated with Generic Heat "
     "Sources, DOE/NETL-4001/09/091014 (compressor 0.85); turbine 0.90 a common sCO2 assumption"),
    ("R11", "eps_recup", "0.95",
     "Effect of actual recuperators' effectiveness on the attainable efficiency of sCO2 Brayton "
     "cycles for solar thermal power plants. Energies 15 (2022) 7773. doi:10.3390/en15207773"),
    ("R12", "rad_sides, rad_areal_mass, T_hab_rad", "2-sided, 14 kg/m2, 283 K",
     "Schunk, R.G. et al. (2022). Thermal Control System Architecture and Technology Challenges "
     "for a Lunar Surface Habitat. IEEE Aerospace Conf., NTRS 20210026557; Schunk et al., "
     "TFAWS 2021, NTRS 20210020015"),
    ("R13", "V_hab, n_crew", "419 m3, 4 crew",
     "Kennedy, K.J. et al. (2008). Lunar Architecture Team Phase 2 Habitat Volume Estimation. "
     "Earth & Space 2008, ASCE. NTRS 20080013517"),
    ("R14", "eps_MLI", "0.02",
     "MLI effective emissivity typically 0.015-0.03: Applied Sciences 12 (2022) 2858. "
     "doi:10.3390/app12062858"),
    ("R15", "f_carnot_hp, T_hp_cond", "0.43, 360 K",
     "Goldman, J., Lovell, T., Ewert, M., Petete, P. (1994). High Temperature Lift Heat Pump "
     "Refrigerant and Thermodynamic Cycle Selection. SAE 941272"),
    ("R16", "q_metab", "14 MJ/person/day",
     "A More Completely Defined CELSS. SAE Technical Paper 941292 (1994)"),
    ("R17", "eta_gen, T_comp_in, pressures", "0.95, 320 K, 8 / 25 MPa", "Paper Table 2"),
    ("R18", "PHX / gas-to-sCO2 heat transfer (context for dT_PHX)", "U ~ 300-400 W/m2-K",
     "Albrecht, K.J., Ho, C.K. et al. (Sandia), Development and testing of a 20 kWt moving "
     "packed-bed particle-to-sCO2 heat exchanger, OSTI 1883488"),
    ("R19", "w_helio_drive", "0.5 W/m2",
     "Heliostat tracking drives: commercial fields ~0.5-4.5 % of plant output (US Patent "
     "9,170,032); 0.5 W/m2 assumed for small lunar heliostats"),
    ("R20", "illumination check (Step 8)", "92.3 % at 2 m, darkness typically 3-5 d",
     "Glaser, P., Scholten, F., De Rosa, D., Marco Figuera, R., Oberst, J., Mazarico, E., Neumann, "
     "G.A., Robinson, M.S. (2014). Illumination conditions at the lunar south pole using high "
     "resolution Digital Terrain Models from LOLA. Icarus 243, 78-90. doi:10.1016/j.icarus.2014.08.013"),
    ("R21", "latitude_deg, horizon below level", "89.47 S (site CR1)",
     "De Rosa, D. et al. (2012). Characterisation of potential landing sites for the European "
     "Space Agency's Lunar Lander project. Planetary and Space Science 74. arXiv:1208.5587. "
     "Horizon method: Mazarico, E. et al. (2011). Icarus 211, 1066-1081. doi:10.1016/j.icarus.2010.10.030"),
    ("R22", "T_surf_night, T_ground", "53 K",
     "Fincannon, H.J., NASA GRC. Lunar Environment and Lunar Power Needs. NTRS 20205002224 "
     "(south-pole candidate sites +22 C to -220 C; sun near the horizon, intermittent darkness)"),
    ("R23", "rho_TES, brick cp check", "2400 kg/m3, ~0.74 J/g-K",
     "Fateri, M., Sottong, R., Kolbe, M., Gamer, J., Sperl, M., Cowley, A. (2019). Thermal "
     "properties of processed lunar regolith simulant. Int. J. Applied Ceramic Technology. "
     "doi:10.1111/ijac.13267 (sintered simulant 2300-2570 kg/m3, k 0.55-1.2 W/m-K)"),
]

ASSUMED = [
    ("Two TES gas loops in counterflow (paper Fig. 4)", "charging: receiver -> end A -> end B (sun only); discharging: end B -> end A -> PHX"),
    ("P_loop = 15 MPa (cycle stays at 25 MPa)", "TES gas loops, sCO2; brick containment at 15 MPa to be justified"),
    ("T_rec_set = 810 K, held by blower 1", "receiver gas outlet; blower 1 speed varies with the sun"),
    ("dT_charge = 150 K", "smallest gas rise in the receiver: caps blower 1 flow when end B is hot"),
    ("N_segments = 10, ideal exchange per segment", "gas leaves each segment at its brick temperature; no conduction between segments"),
    ("dp_seg = 0.25 bar, dp_other_frac = 1 %", "brick and piping pressure drop per loop"),
    ("eta_blower = 0.70, dT_wall = 20 K", "both blowers; receiver wall above the mean gas temperature"),
    ("dT_PHX = 10 K", "PHX (gas -> cycle sCO2): hot-end approach and smallest gap"),
    ("Coolant: water, 1.5 MPa, <= 450 K", "coolant loop fluid and limits"),
    ("dT_precool = 10 K, dp_cool_loop = 150 kPa, eta_pump = 0.7", "precooler gap and coolant pump"),
    ("f_intercept = 0.95, beta_target_deg = 30", "heliostat spillage/blocking, field geometry"),
    ("dp_PHX, dp_recup, dp_cooler = 1 %", "component pressure drops (cycle side)"),
    ("Recuperator always in the flow path", "heat flows backwards when the turbine exhaust is colder than the compressor outlet"),
    ("W_aux_fixed = 3 kWe", "controls, sensors, power electronics"),
    ("TES_mass_margin = 1.20, rad_area_margin = 1.00", "design margins; radiator = peak need of normal operation"),
    ("L_ins = 1 m", "loose-regolith insulation around the TES"),
    ("D_hab = 6 m, f_elec_inside = 1.0", "habitat shape; all 100 kWe dissipated inside"),
    ("Radiator sees deep space only", "sink 3 K at all times; fin efficiency 1"),
    ("Synthetic horizon (dip 2.0 deg, one 4.0 deg range)", "calibrated to south-pole ridge illumination [R20, R21]; replace with a LOLA horizon"),
    ("Design month = darkest synodic month, repeated", "sizing case: the longest darkness of the year every month (conservative)"),
    ("heliostat_spacing_ratio = None", "no mutual shading of heliostats at low sun (optimistic at the pole)"),
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
# ONE results PDF: summary page, report, figures, references (run log only if PDF_RUN_LOG)
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
        fig.text(0.06, 0.965, title + ('  (continued)' if k else ''), fontsize=11, weight='bold')
        fig.text(0.06, 0.945, '\n'.join(wrapped[k:k + per_page]), fontsize=fontsize,
                 family='monospace', va='top')
        pdf.savefig(fig)
        plt.close(fig)
        n += 1
    return n

def _summary_lines():
    s = [f"Run {datetime.datetime.now():%Y-%m-%d %H:%M} | lunar south pole, latitude "
         f"{latitude_deg:.2f} deg | time step {dt:g} h",
         f"Sun: lit {100*SITE_STATS['illum_fraction']:.1f} % of the year | design month {t_lunar:.0f} h: "
         f"sunlit {t_day:.0f} h, longest darkness {t_night:.0f} h",
         f"Demand {W_elec_target:.0f} kWe | TES {T_TES_max:.0f} -> {T_TES_min:.0f} K, {N_segments} segments, "
         f"gas loops {P_loop/1e6:.0f} MPa | cycle {P_low/1e6:.0f} -> {P_high/1e6:.0f} MPa",
         "", f"{'':40}" + "".join(f"{ct.capitalize():>16}" for ct in cycles)]
    def row(name, fn, fmt, unit=''):
        s.append(f"{name:<40}" + "".join(f"{format(fn(ct), fmt):>16}" for ct in cycles) + f"  {unit}")
    row("Cycle efficiency at full charge", lambda ct: cycle_at(TIT_max, ct)['eta'] * 100, '.1f', '%')
    row("TES mass", lambda ct: R[ct]['m_TES'] / 1e3, ',.0f', 't')
    row("TES volume", lambda ct: R[ct]['V_TES'], ',.0f', 'm3')
    row("Collector (mirror) area", lambda ct: R[ct]['A_col'], ',.0f', 'm2')
    row("Radiator panel area", lambda ct: A_RAD[ct], ',.0f', 'm2')
    row("Radiator mass", lambda ct: A_RAD[ct] * rad_areal_mass / 1e3, '.1f', 't')
    row("Gas from end A, lowest", lambda ct: R[ct]['cycle']['T_A'].min(), '.0f', 'K')
    row("Solar heat curtailed", lambda ct: _pct(_tot(R[ct]['cycle'], 'Q_curtail'),
        _tot(R[ct]['cycle'], 'Q_solar')), '.1f', '%')
    row("Parasitic share of gross output", lambda ct: _pct(_tot(R[ct]['cycle'], 'W_aux'),
        _tot(R[ct]['cycle'], 'W_gross')), '.1f', '%')
    row("Demand served (energy)", lambda ct: _pct(_tot(R[ct]['cycle'], 'W_elec'),
        W_elec_target * t_lunar / 1e3), '.1f', '%')
    row("Light on receiver -> electricity", lambda ct: _pct(_tot(R[ct]['cycle'], 'W_elec'),
        _tot(R[ct]['cycle'], 'Q_sun_in')), '.1f', '%')
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
        img = plt.imread(io.BytesIO(FIG_PNG[k - 1]))
        page = plt.figure(figsize=(8.5, 11))
        w_in = 2 * FIG_W
        h_in = w_in * img.shape[0] / img.shape[1]
        ax_img = page.add_axes([(8.5 - w_in) / 2 / 8.5, 1 - (1.0 + h_in) / 11, w_in / 8.5, h_in / 11])
        ax_img.imshow(img, interpolation='none')
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