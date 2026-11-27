# lunar_solar.py -- sun at the site, light on the receiver, receiver loss, surroundings
# Loaded by Step 8 of the notebook (keep this file next to the notebook).
# Not a stand-alone script: it uses the inputs (Step 2) and np, say() from Step 1.
# =====================================================================
# 1. SUN POSITION (standard spherical astronomy, e.g. Meeus, Astronomical Algorithms)
#    Hour angle  w   = 2 pi t / t_lunar - pi/2        (one turn per synodic month)
#    Declination dec = obliquity * sin(2 pi t / t_draconic + season_phase)   (polar seasons)
#    Unit vector to the sun (east, north, up) at latitude phi:
#      e = -cos(dec) sin(w),  n = cos(phi) sin(dec) - sin(phi) cos(dec) cos(w),
#      u =  sin(phi) sin(dec) + cos(phi) cos(dec) cos(w) = sin(elevation)
# 2. HORIZON (horizon method: Mazarico et al. 2011; Glaser et al. 2014)
#    Local horizon elevation vs azimuth from the table in Step 2. Visible fraction of the
#    solar disk (radius 0.2665 deg), x = (elevation - horizon) / radius, clipped to [-1, 1]:
#      visible = 1 - (arccos(x) - x sqrt(1 - x^2)) / pi
# 3. HELIOSTATS: mirror normal bisects the sun and receiver directions,
#    cosine factor = cos(0.5 arccos(s.r)), r at beta_target_deg elevation; averaged over
#    n_field_az mirrors all round the tower (the sun circles the horizon at the pole).
#    Optional low-sun shading: f_shade = min(1, D/H tan(elevation)) (None = no shading).
# 4. RECEIVER
#    light on receiver = G_sun * visible * cosine factor * f_shade * eta_opt   (per m2 mirror)
#    to the gas        = alpha_rec * light - eps_rec sigma (T_rec^4 - T_space^4) / C_ratio, >= 0
# 5. SURROUNDINGS
#    ground (habitat walls): T_surface = max(T_surf_night, T_surf_noon (visible sin(el))^0.25)
#    radiator sink: T_space (vertical panels shaded from the ground see deep space only)
# 6. DESIGN MONTH: the simulation repeats the darkest synodic month of the year, with the
#    longest darkness in the middle; its small seasonal drift is removed so it repeats smoothly.
#    A time step is DARK when the sun is hidden at its midpoint.
# =====================================================================
SUN_RADIUS_DEG = 0.2665            # deg, angular radius of the sun seen from the Moon

# --- 0. Simulation time -> site time on the design month (set in section 6) -------------
T_WIN0, P_WIN, _DEC_DRIFT = 0.0, n_cycle * dt, 0.0
_ABS = False                       # True only for year statistics (no mapping)

def site_time(t):
    t = np.asarray(t, dtype=float)
    return t if _ABS else T_WIN0 + np.mod(t, P_WIN)

# --- 1-2. Sun position and visible fraction ---------------------------------------------
def _sun_site(ts, detrend):
    """Elevation (deg), azimuth (deg), visible fraction at site time ts."""
    ts = np.asarray(ts, dtype=float)
    w = 2 * np.pi * ts / t_lunar - np.pi / 2
    dec = np.radians(obliquity_deg) * np.sin(2 * np.pi * ts / t_draconic + np.radians(season_phase_deg))
    if detrend:
        dec = dec - _DEC_DRIFT * (ts - T_WIN0) / P_WIN
    phi = np.radians(latitude_deg)
    e = -np.cos(dec) * np.sin(w)
    n = np.cos(phi) * np.sin(dec) - np.sin(phi) * np.cos(dec) * np.cos(w)
    u = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(w)
    el = np.degrees(np.arcsin(np.clip(u, -1, 1)))
    az = np.mod(np.degrees(np.arctan2(e, n)), 360.0)
    d = el - horizon_elevation(az)
    x = np.clip(d / SUN_RADIUS_DEG, -1, 1)
    vis = 1 - (np.arccos(x) - x * np.sqrt(1 - x * x)) / np.pi
    return el, az, vis, (e, n, u)

def horizon_elevation(az_deg):
    """Local horizon elevation (deg) at azimuth az_deg."""
    return np.interp(np.mod(az_deg, 360.0), horizon_az_deg, horizon_el_deg, period=360.0)

def sun_state(t):
    """Sun elevation (deg), azimuth (deg) and visible fraction (0-1) at simulation time t."""
    el, az, vis, _ = _sun_site(site_time(t), detrend=not _ABS)
    return el, az, vis

def sun_visible(t):
    return sun_state(t)[2]

# --- 3. Heliostat field -----------------------------------------------------------------
def cosine_factor(t):
    """Heliostat cosine factor, averaged over the field (0 when the sun is hidden)."""
    _, _, vis, (e, n, u) = _sun_site(site_time(t), detrend=not _ABS)
    b = np.radians(beta_target_deg)
    a = np.linspace(0, 2 * np.pi, int(n_field_az), endpoint=False)    # mirror directions
    e, n, u = (np.asarray(x, dtype=float)[..., None] for x in (e, n, u))
    sr = np.cos(b) * (e * np.sin(a) + n * np.cos(a)) + u * np.sin(b)
    cf = np.cos(0.5 * np.arccos(np.clip(sr, -1, 1))).mean(axis=-1)
    return np.where(vis > 0, cf, 0.0)

def shading_factor(t):
    if heliostat_spacing_ratio is None:
        return np.ones_like(np.asarray(t, dtype=float))
    el = sun_state(t)[0]
    return np.minimum(1.0, heliostat_spacing_ratio * np.tan(np.radians(np.maximum(el, 0.0))))

# --- 4. Receiver ------------------------------------------------------------------------
def solar_intercepted(t):
    """Concentrated sunlight reaching the receiver, per m2 of mirror, W/m2."""
    vis = sun_visible(t)
    return G_sun * eta_opt * cosine_factor(t) * shading_factor(t) * vis

def receiver_loss(T_rec):
    """Radiation loss from the receiver aperture, per m2 of mirror, W/m2."""
    return eps_rec * sigma * (np.asarray(T_rec, dtype=float)**4 - T_space**4) / C_ratio

def solar_flux(t, T_rec):
    """Net heat into the charging gas, per m2 of mirror, W/m2 (0 when the receiver would lose)."""
    absorbed = alpha_rec * solar_intercepted(t)
    return np.where(absorbed > 0, np.maximum(absorbed - receiver_loss(T_rec), 0.0), 0.0)

# --- 5. Surroundings --------------------------------------------------------------------
def T_surface(t):
    """Flat ground temperature (K), in balance with the visible sun."""
    el, _, vis = sun_state(t)
    s = np.where(vis > 0, np.maximum(np.sin(np.radians(el)), 0.0), 0.0) * vis
    return np.maximum(T_surf_night, T_surf_noon * s**0.25)

def T_sink_t(t):
    """Radiator sink temperature (K): deep space."""
    return np.full_like(np.asarray(t, dtype=float), T_sink_night)

# --- Speed: each function is calculated once per time step and reused -------------------
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
    wrapper.clear = store.clear
    return wrapper

sun_state, sun_visible, cosine_factor, solar_intercepted, T_surface, T_sink_t = (
    _once_per_time(f) for f in (sun_state, sun_visible, cosine_factor, solar_intercepted,
                                 T_surface, T_sink_t))

# --- 6. Design month --------------------------------------------------------------------
def _runs(mask):
    """(start, length) of each run of True."""
    m = np.diff(np.concatenate(([0], np.asarray(mask, int), [0])))
    st, en = np.where(m == 1)[0], np.where(m == -1)[0]
    return st, en - st

# longest darkness over one draconic year (0.125 h grid) -> month centred on it
_g = 0.125
_tt = np.arange(0.0, t_draconic + 2 * P_WIN, _g)
_st, _ln = _runs(_sun_site(_tt, detrend=False)[2] <= 0)
if _ln.size == 0:
    raise RuntimeError("The site is never dark: nothing for the TES to carry")
T_WIN0 = _st[np.argmax(_ln)] * _g - (n_cycle // 2) * dt
_DEC_DRIFT = float(np.radians(obliquity_deg) * (
    np.sin(2 * np.pi * (T_WIN0 + P_WIN) / t_draconic + np.radians(season_phase_deg))
    - np.sin(2 * np.pi * T_WIN0 / t_draconic + np.radians(season_phase_deg))))
# dark steps of the month (sun hidden at mid-step) and the longest dark run
_dark_steps = _sun_site(T_WIN0 + (np.arange(n_cycle) + 0.5) * dt, detrend=True)[2] <= 0
_st, _ln = _runs(np.concatenate((_dark_steps, _dark_steps)))     # twice: runs across the wrap
_keep = _st < n_cycle
_j = np.argmax(np.minimum(_ln[_keep], n_cycle))
i_dark0, n_dark = int(_st[_keep][_j]), int(min(_ln[_keep][_j], n_cycle))
n_lit = n_cycle - n_dark
t_dark0 = i_dark0 * dt                          # h, start of the longest darkness
t_rise = ((i_dark0 + n_dark) % n_cycle) * dt    # h, its end ("sunrise")
n_night, n_day = n_dark, n_lit
t_night, t_day = n_dark * dt, n_lit * dt
for _f in (sun_state, sun_visible, cosine_factor, solar_intercepted, T_surface, T_sink_t):
    _f.clear()

def is_dark_step(t):
    """True if the sun is hidden at the middle of the step that starts at t."""
    return bool(sun_visible(t + 0.5 * dt) <= 0)

# --- Site statistics over one draconic year (for the report) -----------------------------
_ABS = True
_yr = np.arange(0.0, t_draconic, 1.0)
_el, _, _vis = sun_state(_yr)
_cf = cosine_factor(_yr)
_ABS = False
_lit = _vis > 0
SITE_STATS = dict(illum_fraction=float(_vis.mean()),
                  longest_dark_h=float(_runs(~_lit)[1].max() if (~_lit).any() else 0.0),
                  el_min_deg=float(_el.min()), el_max_deg=float(_el.max()),
                  mean_cosine_lit=float(_cf[_lit].mean()) if _lit.any() else float('nan'),
                  T_WIN0=T_WIN0, t_dark0=t_dark0, t_rise=t_rise, n_dark=n_dark,
                  dark_hours_month=float(_dark_steps.sum() * dt),
                  lit_fraction_month=float(1 - _dark_steps.mean()))

say(f"South pole: sun lit {100*SITE_STATS['illum_fraction']:.1f} % of the year, longest darkness "
    f"{SITE_STATS['longest_dark_h']:.0f} h, sun elevation {SITE_STATS['el_min_deg']:.1f} to "
    f"{SITE_STATS['el_max_deg']:.1f} deg, mean cosine factor {SITE_STATS['mean_cosine_lit']:.3f}")
say(f"Design month: {P_WIN:.0f} h, longest darkness {t_night:.0f} h (t = {t_dark0:.0f}-"
    f"{t_dark0 + t_night:.0f} h), sunrise t = {t_rise:.0f} h, dark in total "
    f"{SITE_STATS['dark_hours_month']:.0f} h | radiator sink {T_sink_night:.0f} K")