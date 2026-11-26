# lunar_solar.py -- sun position, solar input, receiver losses and radiator environment
# Loaded by Step 8 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses the inputs (Step 2) and np, say() from Step 1.
#
# Version 2 (2026-11-26): general sun geometry for any site, including the lunar south pole.
# With the defaults below (latitude 0, obliquity 0, flat horizon, point sun, pole-side
# heliostat, vertical radiator) it reproduces version 1 exactly.
# Speed-up (2026-11-26): the sun is calculated once per time and reused (same numbers).
#
# =====================================================================
# 1. SUN POSITION
# ---------------------------------------------------------------------
# t = hours since the start of the run. Hour angle w advances one full turn per synodic
# month and is -90 deg at t = 0 (sunrise on flat ground at the equator):
#     w = 2 pi t / t_lunar - pi/2
# Sun declination (subsolar latitude) swings with the lunar obliquity over the draconic
# year (the polar "seasons"; Glaser et al. 2017: at the poles the sun moves only
# +-1.54 deg about the horizon):
#     dec = obliquity * sin(2 pi t / t_draconic + season_phase)
# Unit vector to the sun (east, north, up) at site latitude phi (standard spherical
# astronomy, e.g. Meeus, Astronomical Algorithms):
#     e = -cos(dec) sin(w)
#     n =  cos(phi) sin(dec) - sin(phi) cos(dec) cos(w)
#     u =  sin(phi) sin(dec) + cos(phi) cos(dec) cos(w)      = sin(elevation)
#     azimuth = atan2(e, n)   (0 = north, 90 deg = east)
# With dec = 0 this is the version-1 vector ( -sin w, -sin(phi) cos w, cos(phi) cos w ).
#
# 2. TERRAIN SHADOWING (horizon method: Mazarico et al. 2011; Glaser et al. 2014)
# ---------------------------------------------------------------------
# The local horizon is a table of horizon elevation vs azimuth (horizon_az_deg,
# horizon_el_deg). Flat ground = 0 deg everywhere. For a polar site use a published
# profile (Fincannon 2008) or one computed from LOLA topography.
# Visible fraction of the solar disk (angular radius 0.2665 deg), with
# x = (elevation - horizon) / disk radius, clipped to [-1, 1]:
#     visible = 1 - (arccos(x) - x sqrt(1 - x^2)) / pi       (circular segment)
# (as in Raya-Armenta et al. 2023). sun_disk = False treats the sun as a point:
#     visible = 1 if elevation > horizon else 0               (version 1)
#
# 3. HELIOSTAT FIELD
# ---------------------------------------------------------------------
# Mirror normal bisects the sun and receiver directions: cosine factor = cos(0.5 arccos(s.r)).
# r = unit vector from a heliostat to the receiver, at beta_target_deg elevation.
#   field_layout = 'pole_side': one representative heliostat on the pole side of the tower,
#                  r points toward the equator (version 1).
#   field_layout = 'surround' : n_field_az heliostats evenly around the tower, cosine factor
#                  averaged over the field. Use this at the pole, where the sun circles the
#                  whole horizon each month.
# Low-sun shading/blocking between mirrors (shadow of a mirror of height H on the next one
# at ground spacing D):  f_shade = min(1, D/H * tan(elevation)).
#   heliostat_spacing_ratio = None -> no shading (version 1).
#
# 4. RECEIVER (unchanged)
# ---------------------------------------------------------------------
#   absorbed  = G_sun * visible * cosine factor * f_shade * eta_opt * alpha_rec
#   rec. loss = eps_rec * sigma * (T_rec^4 - T_space^4) / C_ratio
#   to gas    = absorbed - loss, never below 0 (receiver closed when it would lose heat)
# T_rec = mean gas temperature in the receiver + dT_wall, from the TES gas loop each hour
# (Steps 6 and 9).
#
# 5. SURFACE AND RADIATOR SINK
# ---------------------------------------------------------------------
# Flat ground in radiative balance with the visible sun (no thermal inertia):
#     T_surface = max(T_surf_night, T_surf_noon * (visible * sin(elevation))^0.25)
# T_surf_noon = 392 K is the equatorial subsolar value, so the same law gives the cold
# polar ground when the sun is low.
# Radiator, effective sink (per radiating face):
#     T_sink^4 = F * T_surface^4 + (1 - F) * T_space^4 + q_sun / (emissivity * sigma)
#   radiator_orientation = 'vertical'  : F = F_rad_surf, panel edge-on to the sun, q_sun = 0
#                                        (version 1)
#   radiator_orientation = 'horizontal': upward face sees only sky (F = 0) and absorbs
#                                        q_sun = alpha_rad * G_sun * visible * sin(elevation)
#   shield = True: F = 0 for either orientation.
# Not included: sunlight reflected by the surface, regolith thermal inertia, local slopes.
# =====================================================================

# --- New inputs (set them in Step 2; defaults reproduce version 1) -------------------
_v2_defaults = {
    'obliquity_deg': 0.0,            # deg, lunar obliquity; 1.54 for real seasons
    't_draconic': 346.62 * 24,       # h, draconic year (period of the declination swing)
    'season_phase_deg': 0.0,         # deg, phase of the declination at t = 0
                                     #   90 -> dec = +obliquity at t = 0 (southern winter, worst
                                     #   case for a south-pole site); -90 -> northern winter
    'horizon_az_deg': np.array([0.0, 360.0]),   # deg, azimuths of the horizon table
    'horizon_el_deg': np.array([0.0, 0.0]),     # deg, horizon elevation (flat ground)
    'sun_disk': False,               # True: partial sun at the horizon (finite solar disk)
    'field_layout': 'pole_side',     # 'pole_side' (version 1) or 'surround'
    'n_field_az': 36,                # heliostat directions around the tower ('surround')
    'heliostat_spacing_ratio': None, # D/H, ground spacing / mirror height; None = no shading
    'radiator_orientation': 'vertical',   # 'vertical' (version 1) or 'horizontal'
    'alpha_rad': 0.20,               # radiator solar absorptance ('horizontal' only)
}
for _k, _v in _v2_defaults.items():
    if _k not in globals():
        globals()[_k] = _v

SUN_RADIUS_DEG = 0.2665             # deg, angular radius of the sun seen from the Moon
_EL_TOL_DEG = 1e-9                  # deg, round-off guard: sun exactly on the horizon = not up

# --- 1. Sun position ------------------------------------------------------------------
def _hour_angle(t):
    """Hour angle (rad): -pi/2 at t = 0, one full turn per synodic month."""
    return 2.0 * np.pi * np.asarray(t, dtype=float) / t_lunar - np.pi / 2

def sun_declination(t):
    """Sun declination (rad) at the Moon: +-obliquity over the draconic year."""
    t = np.asarray(t, dtype=float)
    return np.radians(obliquity_deg) * np.sin(2.0 * np.pi * t / t_draconic
                                              + np.radians(season_phase_deg))

def sun_vector(t):
    """Unit vector to the sun (east, north, up) at the site."""
    w, dec, phi = _hour_angle(t), sun_declination(t), np.radians(latitude_deg)
    e = -np.cos(dec) * np.sin(w)
    n = np.cos(phi) * np.sin(dec) - np.sin(phi) * np.cos(dec) * np.cos(w)
    u = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(w)
    return e, n, u

def horizon_elevation(az_deg):
    """Local horizon elevation (deg) at azimuth az_deg (deg), from the horizon table."""
    return np.interp(np.mod(az_deg, 360.0), horizon_az_deg, horizon_el_deg, period=360.0)

def sun_state(t):
    """Sun elevation (deg), azimuth (deg) and visible fraction of the solar disk (0-1)."""
    e, n, u = sun_vector(t)
    el = np.degrees(np.arcsin(np.clip(u, -1.0, 1.0)))
    az = np.mod(np.degrees(np.arctan2(e, n)), 360.0)
    d = el - horizon_elevation(az)
    if sun_disk:
        x = np.clip(d / SUN_RADIUS_DEG, -1.0, 1.0)
        vis = 1.0 - (np.arccos(x) - x * np.sqrt(1.0 - x * x)) / np.pi
    else:
        vis = (d > _EL_TOL_DEG).astype(float)
    return el, az, vis

def sun_sin_elevation(t):
    """sin(sun elevation) while the sun is visible, else 0 (name kept from version 1)."""
    el, _, vis = sun_state(t)
    return np.where(vis > 0, np.maximum(np.sin(np.radians(el)), 0.0), 0.0)

def sun_visible(t):
    """Visible fraction of the solar disk (0-1). Use this, not t < t_day, to tell day."""
    return sun_state(t)[2]

# --- 3. Heliostat field ---------------------------------------------------------------
def _cos_half_angle(sr):
    return np.cos(0.5 * np.arccos(np.clip(sr, -1.0, 1.0)))

def cosine_factor(t):
    """Heliostat cosine factor, averaged over the field (0 when the sun is not visible)."""
    e, n, u = sun_vector(t)
    _, _, vis = sun_state(t)
    b = np.radians(beta_target_deg)
    if field_layout == 'pole_side':
        a_r = 0.0 if latitude_deg <= 0 else np.pi          # r toward the equator
        sr = np.cos(b) * (e * np.sin(a_r) + n * np.cos(a_r)) + u * np.sin(b)
        cf = _cos_half_angle(sr)
    elif field_layout == 'surround':
        a_r = np.linspace(0.0, 2.0 * np.pi, int(n_field_az), endpoint=False)
        e_, n_, u_ = (np.asarray(x, dtype=float)[..., None] for x in (e, n, u))
        sr = np.cos(b) * (e_ * np.sin(a_r) + n_ * np.cos(a_r)) + u_ * np.sin(b)
        cf = _cos_half_angle(sr).mean(axis=-1)
    else:
        raise ValueError("field_layout must be 'pole_side' or 'surround'")
    return np.where(vis > 0, cf, 0.0)

def shading_factor(t):
    """Fraction of mirror area not shaded by neighbours at low sun (1 = no shading)."""
    if heliostat_spacing_ratio is None:
        return np.ones_like(np.asarray(t, dtype=float))
    el, _, _ = sun_state(t)
    tan_el = np.tan(np.radians(np.maximum(el, 0.0)))
    return np.minimum(1.0, heliostat_spacing_ratio * tan_el)

# --- 4. Receiver ----------------------------------------------------------------------
def solar_intercepted(t):
    """Concentrated sunlight reaching the receiver, per m2 of collector, W/m2."""
    _, _, vis = sun_state(t)
    s = sun_sin_elevation(t)
    if collector_mode == 'heliostat':
        factor = cosine_factor(t) * shading_factor(t)
    elif collector_mode == 'fixed':
        factor = s
    else:                                                  # 'tracking'
        factor = (s > 0).astype(float)
    if sun_disk:
        factor = factor * vis
    return G_sun * factor * eta_opt

def receiver_loss(T_rec):
    """Radiation loss from the receiver aperture at receiver temperature T_rec, per m2 of
    collector, W/m2."""
    T_r = np.asarray(T_rec, dtype=float)
    return eps_rec * sigma * (T_r**4 - T_space**4) / C_ratio

def solar_flux(t, T_rec):
    """Net heat delivered to the loop gas per m2 of collector, W/m2."""
    absorbed = alpha_rec * solar_intercepted(t)
    return np.where(absorbed > 0, np.maximum(absorbed - receiver_loss(T_rec), 0.0), 0.0)

def receiver_efficiency(t, T_rec):
    """Net heat to the gas / concentrated light reaching the receiver (0 when closed)."""
    q_in = solar_intercepted(t)
    return np.where(q_in > 0, solar_flux(t, T_rec) / np.where(q_in > 0, q_in, 1.0), 0.0)

# --- 5. Surface and radiator sink -----------------------------------------------------
def T_surface(t):
    """Flat lunar ground temperature (K), radiative balance with the visible sun."""
    _, _, vis = sun_state(t)
    s = sun_sin_elevation(t) * (vis if sun_disk else 1.0)
    return np.maximum(T_surf_night, T_surf_noon * s**0.25)

def T_sink_t(t):
    """Effective radiator sink temperature (K)."""
    if radiator_orientation == 'vertical':
        F = 0.0 if shield else F_rad_surf
        q_sun = 0.0
    elif radiator_orientation == 'horizontal':
        F = 0.0
        _, _, vis = sun_state(t)
        q_sun = alpha_rad * G_sun * sun_sin_elevation(t) * (vis if sun_disk else 1.0)
    else:
        raise ValueError("radiator_orientation must be 'vertical' or 'horizontal'")
    return (F * T_surface(t)**4 + (1 - F) * T_space**4 + q_sun / (emissivity * sigma)) ** 0.25

# --- Speed: the sun is calculated once per time and reused -----------------------------
# The hourly simulation asks for the sun many times per hour (light on the receiver on
# every pass of the hourly iteration, surface and radiator-sink temperatures in every
# plant calculation). For a single time t, each function below is calculated once and the
# answer is reused -- the same calculation, so the same numbers. Arrays of times (plots,
# site statistics) are always calculated. Changing an input means re-running from Step 2
# (Step 8 then starts with empty stores).
def _once_per_time(fn):
    store = {}
    def wrapper(t, *args):
        if args or np.ndim(t) != 0:
            return fn(t, *args)
        key = float(t)
        if key not in store:
            if len(store) > 50000:
                store.clear()
            store[key] = fn(t)
        return store[key]
    wrapper.__name__, wrapper.__doc__ = fn.__name__, fn.__doc__
    wrapper.clear = store.clear
    return wrapper

sun_vector = _once_per_time(sun_vector)
sun_state = _once_per_time(sun_state)
sun_sin_elevation = _once_per_time(sun_sin_elevation)
sun_visible = _once_per_time(sun_visible)
cosine_factor = _once_per_time(cosine_factor)
shading_factor = _once_per_time(shading_factor)
solar_intercepted = _once_per_time(solar_intercepted)
T_surface = _once_per_time(T_surface)
T_sink_t = _once_per_time(T_sink_t)

# --- 6. Site statistics (compare with published illumination studies) -----------------
def _longest_run(mask, dt_h):
    """Longest continuous run of True in a boolean series, in hours."""
    best = run = 0
    for m in mask:
        run = run + 1 if m else 0
        best = max(best, run)
    return best * dt_h

def illumination_stats(t_end=None, dt_h=1.0):
    """Illumination fraction, longest darkness/light (h), sun elevation range, mean
    cosine factor and radiator-sink range over [0, t_end] (default: one draconic year)."""
    t_end = t_draconic if t_end is None else t_end
    tt = np.arange(0.0, t_end, dt_h)
    el, _, vis = sun_state(tt)
    lit = vis > 0
    cf = cosine_factor(tt)[lit]
    Ts = T_sink_t(tt)
    return {
        'illum_fraction': float(vis.mean()),
        'longest_dark_h': _longest_run(~lit, dt_h),
        'longest_light_h': _longest_run(lit, dt_h),
        'el_min_deg': float(el.min()), 'el_max_deg': float(el.max()),
        'mean_cosine_lit': float(cf.mean()) if cf.size else float('nan'),
        'T_sink_min': float(Ts.min()), 'T_sink_max': float(Ts.max()),
    }

# --- Quick look -----------------------------------------------------------------------
_st = illumination_stats(t_end=max(t_lunar, t_draconic if obliquity_deg > 0 else t_lunar))
say(f"Site: latitude {latitude_deg:.2f} deg, obliquity {obliquity_deg:.2f} deg, "
    f"horizon {'flat' if not np.any(horizon_el_deg) else 'from table'}, "
    f"field '{field_layout}', radiator '{radiator_orientation}'")
say(f"  Sun lit {100*_st['illum_fraction']:.1f} % of the time | longest darkness "
    f"{_st['longest_dark_h']:.0f} h | longest light {_st['longest_light_h']:.0f} h | "
    f"sun elevation {_st['el_min_deg']:.1f} to {_st['el_max_deg']:.1f} deg")
say(f"  Radiator sink {_st['T_sink_min']:.0f} to {_st['T_sink_max']:.0f} K | mean cosine "
    f"factor while lit {_st['mean_cosine_lit']:.3f}")
say(f"Receiver: alpha {alpha_rec}, emissivity {eps_rec}, concentration {C_ratio:.0f}x, "
    f"T_rec from the gas loop | optics {eta_opt:.3f}")
_T_rec_show = (600.0, 750.0, 900.0)                     # K, receiver temperatures to tabulate
_tt = np.arange(0.0, max(t_lunar, t_draconic if obliquity_deg > 0 else t_lunar), 0.25)
_el_series, _, _vis_series = sun_state(_tt)
say(f"{'sun elev':>9}" + "".join(f"{f'T_rec {T:.0f} K':>14}" for T in _T_rec_show))
for el in (90, 45, 30, 15, 10, 5, 2, 1):
    hits = np.where((np.abs(_el_series - el) < 0.25) & (_vis_series > 0))[0]
    if hits.size == 0:
        continue                                       # the sun never reaches this elevation here
    t_el = _tt[hits[0]]
    say(f"{el:>8}°" + "".join(f"{float(receiver_efficiency(t_el, T))*100:>13.1f}%"
                              for T in _T_rec_show))