# lunar_regolith.py -- regolith TES properties of the lunar regolith-TES sCO2 model
# Loaded by Step 3 of Lunar_Regolith_TES_Cycle.ipynb (keep this file next to the notebook).
# Not a stand-alone script: it uses np and the inputs (regolith_type, T_TES_min, T_TES_max,
# rho_TES, T_ground, L_ins) set by Step 2.
# ---------------------------------------------------------------------
# Regolith specific heat (Table 3, Schreiner et al. 2016)
# For T > 350 K: cp = Cf + Cg*T + Ch/T^2   [J/kg-K]
# (The 100-350 K polynomial in Table 3 gives J/g-K, not J/kg-K;
#  it is not needed here because the TES stays at or above 350 K.)
# ---------------------------------------------------------------------
REGOLITH_CP = {
    'high_Ti_mare': (8.820e2, 3.083e-1, -2.278e7),
    'low_Ti_mare':  (9.093e2, 2.870e-1, -2.469e7),
    'highlands':    (9.530e2, 2.524e-1, -2.645e7),
}
Cf, Cg, Ch = REGOLITH_CP[regolith_type]

def cp_regolith(T):
    return Cf + Cg * T - (-Ch) / T**2      # J/kg-K  (Ch is negative)

def e_regolith(T):
    """Specific sensible energy relative to T_TES_min, J/kg (integral of cp dT)."""
    F = lambda T: Cf * T + 0.5 * Cg * T**2 - Ch / T
    return F(T) - F(T_TES_min)

# Lookup table for inverting energy -> temperature
T_grid_TES = np.linspace(T_TES_min - 50, T_TES_max + 1, 4000)
e_grid_TES = e_regolith(T_grid_TES)

def T_from_e(e):
    return np.interp(e, e_grid_TES, T_grid_TES)

def k_regolith(T):
    """Thermal conductivity of granular regolith, W/m-K (Table 3)."""
    return 0.001561 + 5.426e-11 * np.asarray(T)**3

def TES_heat_loss(T_TES, m_TES):
    """Heat lost through the insulation, kW. TES = cube of volume m/rho."""
    A = 6 * (m_TES / rho_TES) ** (2 / 3)                     # m2
    k = k_regolith(0.5 * (T_TES + T_ground))                 # insulation at its mean temperature
    return k * A / L_ins * (T_TES - T_ground) / 1e3
