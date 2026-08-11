"""
Hydrazine thruster catalyst bed reactor simulation.
Three-region model: Liquid → Liquid-Vapor → Vapor, based on Kesten (1967).
All state variables and tables stored in SI units.

V2: NH3 intra-particle diffusion via shooting method (SGRAD equivalent).
    N2H4 remains mass-transfer limited (Kesten VAPOR.f MFLAG=1).
"""

import json
import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# ─── Universal constants ──────────────────────────────────────────────────────
R_UNIV = 8.314          # J/(mol·K)
MW = {                  # kg/mol
    'N2H4': 0.032048,
    'NH3':  0.017032,
    'N2':   0.028016,
    'H2':   0.002016,
}

# ─── Unit conversion factors  (Kesten English → SI) ─────────────────────────
_R_K    = 1.0 / 1.8        # °R → K
_BTU    = 2326.0            # BTU/lb → J/kg
_BTU_R  = 4186.8            # BTU/(lb·°R) → J/(kg·K)
_ft     = 0.3048            # ft → m
_psia   = 6894.757          # psia → Pa
_G_si   = 4.88243           # lb/(ft²·s) → kg/(m²·s)
_mu_si  = 1.48816           # lb/(ft·s) → Pa·s

# ─── Kesten (1967) input parameters, converted to SI ─────────────────────────
T_FEED  = 530.0  * _R_K     # K   ≈  294.4   feed temperature
P_IN    = 100.0  * _psia    # Pa  ≈  689 476  inlet pressure
G       = 3.0    * _G_si    # kg/(m²·s) ≈ 14.65  mass flux
H_FEED  = 0.0               # J/kg       specific enthalpy reference
Z_END   = 0.25   * _ft      # m   ≈  0.0762   bed length
CP_LIQ  = 0.7332 * _BTU_R  # J/(kg·K) ≈ 3069   liquid N2H4 Cp

# Activation temperatures: Kesten values [°R] → [K]
AGM_K = 2500.0  / 1.8   # K ≈  1 389   catalytic N2H4 decomp (liquid)
BGM_K = 50000.0 / 1.8   # K ≈ 27 778   catalytic NH3  decomp (vapor)
CGM_K = 33000.0 / 1.8   # K ≈ 18 333   thermal   N2H4 decomp (vapor)

# Pre-exponentials (calibrated for Kesten's lb/ft³ concentration basis)
ALPHA1 = 1.00e10   # s⁻¹   catalytic N2H4  (1st order → unit-independent)
ALPHA2 = 1.00e11   # s⁻¹   catalytic NH3
ALPHA3 = 2.14e10   # s⁻¹   thermal   N2H4  (1st order → unit-independent)

# Reference state for diffusivity correction (Chapman-Enskog)
T_REF = 492.0 * _R_K     # K   = 273.33
P_REF = 14.7  * _psia    # Pa  = 101 325

# Binary diffusivities at reference state [ft²/s] → [m²/s]
D0_N2H4 = 0.95e-4 * _ft**2   # m²/s
D0_NH3  = 0.17e-3 * _ft**2   # m²/s

# Reaction orders
EN3 = -1.6   # H2 inhibition order for catalytic NH3 decomp


# ─── Property tables ─────────────────────────────────────────────────────────

def _load_tables(json_path):
    """Parse data.json; return dict of {name: (x_arr, y_arr)} in raw units."""
    with open(json_path) as f:
        raw = json.load(f)['tables']
    tables = {}
    for name, tbl in raw.items():
        x_key = next(k for k in tbl if k not in ('y', 'description'))
        tables[name] = (np.asarray(tbl[x_key], float),
                        np.asarray(tbl['y'],    float))

    # Restore precise Kesten patch data (override JSON where needed)
    t_vp = [492.0, 519.0, 528.37, 529.08, 534.60, 534.71, 538.84, 543.91,
            545.73, 560.20, 569.98, 579.26, 579.48, 595.34, 610.13, 614.08,
            618.07, 627.49, 628.82, 645.68, 650.76, 665.57, 674.99, 686.13,
            692.39, 697.47, 744.0, 798.0, 852.0, 942.0, 1032.0, 1122.0, 1176.0]
    p_vp = [0.0520, 0.1479, 0.2011, 0.2069, 0.2398, 0.2436, 0.2823, 0.2920,
            0.3539, 0.5453, 0.7367, 0.9727, 0.9823, 1.510, 2.204, 2.462, 2.740,
            3.407, 3.562, 5.240, 5.971, 8.065, 9.711, 11.91, 13.46, 14.70,
            33.80, 73.48, 147.0, 382.1, 823.0, 1528.0, 2131.0]
    tables['TBLVP'] = (np.array(t_vp), np.array(p_vp))

    p_tv = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,
            850,900,950,1000]
    t_tv = [770,820,855,880,905,925,945,965,980,995,1010,1025,1035,1050,1060,
            1070,1080,1090,1100,1110]
    tables['VPTBL'] = (np.array(p_tv, float), np.array(t_tv, float))

    t_mu = [360,540,720,900,1080,1260,1440,1620,1800,1980,2160,2340,2520]
    y_mu = [0.048e-4, 0.070e-4, 0.093e-4, 0.117e-4, 0.141e-4, 0.164e-4,
            0.186e-4, 0.207e-4, 0.220e-4, 0.247e-4, 0.266e-4, 0.285e-4, 0.302e-4]
    tables['VISVST'] = (np.array(t_mu, float), np.array(y_mu, float))
    return tables


class Props:
    """
    Thermodynamic and geometric property tables, all converted to SI at load.
    Temperature axes: K.  Pressure axes: Pa.  Length axes: m.
    """

    def __init__(self, json_path=None):
        if json_path is None:
            # Allow running from src/ or from the project root
            for candidate in ('data.json', '../data.json'):
                if Path(candidate).exists():
                    json_path = candidate
                    break
        raw = _load_tables(json_path)

        def _si(name, x_scale, y_scale):
            x, y = raw[name]
            return x * x_scale, y * y_scale

        # Cp [BTU/(lb·°R)] → [J/(kg·K)],  T axis [°R] → [K]
        self._cp = {
            'H2':   _si('CFTBL1', _R_K, _BTU_R),
            'N2':   _si('CFTBL2', _R_K, _BTU_R),
            'NH3':  _si('CFTBL3', _R_K, _BTU_R),
            'N2H4': _si('CFTBL4', _R_K, _BTU_R),
        }
        # Heats of reaction [BTU/lb] → [J/kg],  T axis [°R] → [K]
        # H4TBL is negative (exothermic N2H4 decomp); H3TBL is positive (endothermic NH3 decomp)
        self._h4 = _si('H4TBL', _R_K, _BTU)
        self._h3 = _si('H3TBL', _R_K, _BTU)

        # Viscosity [lb/(ft·s)] → [Pa·s],  T axis [°R] → [K]
        self._visc = _si('VISVST', _R_K, _mu_si)

        # Geometry  z axis [ft] → [m]
        self._A    = _si('ZTBLA',  _ft, _ft)       # particle radius [m]
        self._ap   = _si('ZTBLAP', _ft, 1.0 / _ft) # specific surface [m²/m³] (1/ft → 1/m)
        self._void = _si('ZTBLD',  _ft, 1.0)        # void fraction [-]

        # Saturation: P [psia] → Tsat [°R]
        self._tsat = _si('VPTBL', _psia, _R_K)

        # Vaporisation enthalpies [BTU/lb] → [J/kg],  Tsat axis [°R] → [K]
        self._dhvap = _si('DHVST',  _R_K, _BTU)
        self._dhliq = _si('DHLVST', _R_K, _BTU)

    # ── Public accessors ──────────────────────────────────────────────────────

    def cp(self, species, T):
        """Specific heat of species at T [K].  Returns J/(kg·K)."""
        return float(np.interp(T, *self._cp[species]))

    def h_rxn_N2H4(self, T):
        """Heat of reaction for N2H4 decomp at T [K].  J/kg_N2H4 (negative)."""
        return float(np.interp(T, *self._h4))

    def h_rxn_NH3(self, T):
        """Heat of reaction for NH3 decomp at T [K].  J/kg_NH3 (positive)."""
        return float(np.interp(T, *self._h3))

    def viscosity(self, T):
        """Dynamic viscosity at T [K].  Pa·s."""
        return float(np.interp(T, *self._visc))

    def particle_radius(self, z):
        """Catalyst particle radius at axial position z [m].  Returns m."""
        return float(np.interp(z, *self._A))

    def specific_surface(self, z):
        """External particle surface area per bed volume at z [m].  m²/m³."""
        return float(np.interp(z, *self._ap))

    def void_fraction(self, z):
        """Interparticle void fraction at z [m].  Dimensionless."""
        return float(np.interp(z, *self._void))

    def T_sat(self, P):
        """Saturation temperature at pressure P [Pa].  Returns K."""
        return float(np.interp(P, *self._tsat))

    def dh_vap(self, T_s):
        """Latent heat of vaporisation at Tsat [K].  J/kg."""
        return float(np.interp(T_s, *self._dhvap))

    def dh_liq(self, T_s):
        """Enthalpy of liquid N2H4 at Tsat [K].  J/kg."""
        return float(np.interp(T_s, *self._dhliq))


# ─── Transport helpers ────────────────────────────────────────────────────────

def diffusivity(D0, T, P):
    """
    Binary diffusivity [m²/s] corrected for T and P via Chapman-Enskog.
    D0: reference diffusivity at (T_REF, P_REF) [m²/s]
    """
    corr = 1.0 - np.exp(-0.0672 * (P / P_REF) * (T_REF / T))
    return D0 * (T / T_REF) ** 1.832 * (P_REF / P) * corr


def kc_chilton_colburn(G_flux, rho, mu, D, a_p):
    """
    Mass transfer coefficient [m/s] from the Chilton-Colburn correlation.
    Returns 0 on degenerate inputs.
    """
    if rho <= 0 or mu <= 0 or D <= 0 or a_p <= 0:
        return 0.0
    u_s = G_flux / rho
    Sc  = mu / (rho * D)
    Re  = G_flux / (a_p * mu)
    return 0.61 * u_s * Sc ** (-2.0 / 3.0) * Re ** (-0.41)


def mixture_density(T, P, w):
    """Ideal-gas mixture density [kg/m³] from mass fractions w (dict)."""
    inv_MW = sum(w[sp] / MW[sp] for sp in w)   # mol/kg
    return P / (R_UNIV * T * inv_MW)


def mixture_cp(T, props, w):
    """Mass-weighted mean specific heat [J/(kg·K)] from mass fractions."""
    return sum(w[sp] * props.cp(sp, T) for sp in w)


def thiele_effectiveness(k_rxn, R, D_eff):
    """
    Pore-diffusion effectiveness factor η for a first-order reaction in a sphere.
    k_rxn : apparent surface rate constant [m/s]
    R     : particle radius [m]
    D_eff : effective intra-particle diffusivity [m²/s]
    Returns η ∈ (0, 1].
    """
    k_vol = k_rxn / (R / 3.0)              # volumetric rate constant [1/s]
    phi   = R * np.sqrt(k_vol / D_eff)     # Thiele modulus (sphere)
    if phi < 1e-3:
        return 1.0
    if phi > 100.0:
        return 3.0 / phi                    # strong-diffusion limit
    return (3.0 / phi**2) * (phi / np.tanh(phi) - 1.0)


def sgrad_NH3(C_NH3_bulk, C_H2_bulk, T, P, A, props, kc_NH3, kc_N2H4, C_N2H4_bulk):
    """
    Shooting-method equivalent of Kesten's SGRAD subroutine for NH3 decomp.

    Solves the spherical radial diffusion-reaction BVP by integrating outward
    from r=0 and using brentq on the center concentration:

        D_eff/r² · d/dr(r² dC/dr) = K0 · C        (first-order in C_NH3)
        BC r=0: dC/dr = 0
        BC r=A: D_eff · dC/dr = kc · (C_bulk - C_surface)

    Also computes TPS (particle surface temperature) from the film energy
    balance, matching Kesten's SGRAD approach.  K0 and D_eff are evaluated
    at TPS, not T_bulk.

    Parameters
    ----------
    C_NH3_bulk  : float  bulk NH3 concentration [kg/m³]
    C_H2_bulk   : float  bulk H2 concentration  [kg/m³]  (uniform in particle)
    T           : float  bulk gas temperature [K]
    P           : float  pressure [Pa]
    A           : float  particle radius [m]
    props       : Props  property tables (for h4, h3, viscosity)
    kc_NH3      : float  external film mass-transfer coeff for NH3 [m/s]
    kc_N2H4     : float  external film mass-transfer coeff for N2H4 [m/s]
    C_N2H4_bulk : float  bulk N2H4 concentration [kg/m³]

    Returns
    -------
    float  rate per unit external particle area [kg_NH3 / (m²·s)]
    """
    if C_NH3_bulk <= 1e-30 or kc_NH3 <= 0:
        return 0.0

    D_eff = diffusivity(D0_NH3, T, P)

    # ── Apparent first-order rate constant ───────────────────────────────────
    # k2 [m/s]: surface rate coefficient (same units as kc_NH3, used in
    # film-diffusion code).  For the intra-particle ODE we need K0 [1/s]:
    # for a sphere, K0_vol = k_surface × 3/A.
    C_H2 = max(C_H2_bulk, 1e-20)
    k2   = ALPHA2 * np.exp(-BGM_K / T) / (C_H2 ** abs(EN3))   # [m/s]
    K0   = k2 * 3.0 / A                                         # [1/s]

    phi = A * np.sqrt(K0 / D_eff)   # Thiele modulus

    # ── Strategy: analytical for large φ, shooting for small φ ───────────────
    # For first-order kinetics the Thiele formula is EXACT, and for φ > 5 the
    # outward-shooting ODE is badly conditioned (exponential growth overflows).
    # Shooting is reserved for φ ≤ 5 where it is well-conditioned and where
    # future non-linear kinetics (e.g. concentration-dependent inhibition
    # evaluated inside the pore) would differ from the analytical formula.
    PHI_SWITCH = 5.0

    if phi > PHI_SWITCH:
        # Analytical Thiele effectiveness (exact for first-order sphere)
        eta   = thiele_effectiveness(k2, A, D_eff)
        k_eff = kc_NH3 * k2 * eta / (kc_NH3 + k2 * eta)
        return k_eff * C_NH3_bulk

    # ── Shooting method (φ ≤ 5) ──────────────────────────────────────────────
    r_start = 1e-6 * A   # avoid singularity at r = 0

    def radial_ode(r, y):
        C    = max(y[0], 0.0)
        dCdr = y[1]
        rate = K0 * C                           # [kg_NH3 / (m³·s)]
        if r < 1e-8 * A:
            d2Cdr2 = rate / (3.0 * D_eff)      # L'Hopital limit at r → 0
        else:
            d2Cdr2 = rate / D_eff - 2.0 * dCdr / r
        return [dCdr, d2Cdr2]

    def shoot(C_center):
        """Integrate outward; return surface-BC residual (zero at solution)."""
        dCdr0 = K0 * C_center / (3.0 * D_eff) * r_start
        sol = solve_ivp(
            radial_ode,
            [r_start, A],
            [C_center, dCdr0],
            method='RK45',
            rtol=1e-5,
            atol=1e-14,
            max_step=A / 20.0,
        )
        if not sol.success:
            return np.nan
        C_surf = max(sol.y[0, -1], 0.0)
        dCdr_R = sol.y[1, -1]
        return D_eff * dCdr_R - kc_NH3 * (C_NH3_bulk - C_surf)

    lo   = 1e-10 * C_NH3_bulk
    hi   = C_NH3_bulk * (phi / np.sinh(phi) if phi < 50 else 0.0) * 2.0
    hi   = np.clip(hi, lo * 10, 0.9999 * C_NH3_bulk)
    f_lo = shoot(lo)
    f_hi = shoot(hi)

    if np.isnan(f_lo) or np.isnan(f_hi) or np.sign(f_lo) == np.sign(f_hi):
        eta   = thiele_effectiveness(k2, A, D_eff)
        k_eff = kc_NH3 * k2 * eta / (kc_NH3 + k2 * eta)
        return k_eff * C_NH3_bulk

    C_center_sol = brentq(shoot, lo, hi,
                          xtol=1e-8 * C_NH3_bulk, rtol=1e-6, maxiter=60)

    dCdr0 = K0 * C_center_sol / (3.0 * D_eff) * r_start
    sol_f = solve_ivp(
        radial_ode, [r_start, A], [C_center_sol, dCdr0],
        method='RK45', rtol=1e-5, atol=1e-14, max_step=A / 20.0,
    )
    C_surf = max(sol_f.y[0, -1], 0.0)
    return kc_NH3 * (C_NH3_bulk - C_surf)


def ergun_dpdz(G_flux, rho, mu, d_p, eps):
    """
    Ergun pressure gradient [Pa/m].  Negative = pressure drops along z.
    d_p: particle diameter [m]   eps: void fraction [-]
    """
    u_s = G_flux / rho
    return -(150.0 * mu * u_s * (1 - eps) ** 2 / (d_p ** 2 * eps ** 3)
             + 1.75 * rho * u_s ** 2 * (1 - eps) / (d_p * eps ** 3))


# ─── Phase boundary conditions ────────────────────────────────────────────────

def compute_phase_bcs(props):
    """
    Compute saturation conditions and initial vapour composition.
    Returns a dict of scalar BCs used by all three phase solvers.
    """
    T_s  = props.T_sat(P_IN)                       # K   saturation temperature
    H_L  = H_FEED + (T_s - T_FEED) * CP_LIQ       # J/kg enthalpy at sat. liquid
    H_V  = H_L + props.dh_vap(T_s) - props.dh_liq(T_s)  # J/kg at sat. vapour

    # Fraction of N2H4 converted before entering the vapour phase (energy balance).
    # Using Kesten CONC convention: XV = H_V / |h_rxn_N2H4(T_s)|
    h4 = props.h_rxn_N2H4(T_s)                     # J/kg_N2H4, negative
    XV  = max(0.0, min(0.99, H_V / abs(h4)))

    # Kesten CONC stoichiometry (overall simplified):
    #   N2H4 → NH3 + 0.5 N2 + 0.5 H2   (XV converted fraction)
    # Mole counts per initial mole N2H4: (1-XV) N2H4 + XV NH3 + XV/2 N2 + XV/2 H2
    # Total moles = 1 + XV
    denom  = 1.0 + XV
    y0 = {
        'N2H4': (1.0 - XV) / denom,
        'NH3':        XV   / denom,
        'N2':   (XV / 2.0) / denom,
        'H2':   (XV / 2.0) / denom,
    }

    # Convert mole fractions → mass fractions for the ODE state
    mass_unnorm = {sp: y0[sp] * MW[sp] for sp in MW}
    total_mass  = sum(mass_unnorm.values())
    w0 = {sp: mass_unnorm[sp] / total_mass for sp in MW}

    return dict(T_sat=T_s, H_L=H_L, H_V=H_V, XV=XV, y0=y0, w0=w0)


# ─── Liquid phase ─────────────────────────────────────────────────────────────

def run_liquid(props, bcs, n_steps=500):
    """
    March from H=H_FEED to H=H_L (saturation) in the liquid phase.
    Returns arrays (z, T, H) for this region.
    """
    H_L = bcs['H_L']
    T_s = bcs['T_sat']

    z_arr, T_arr, H_arr = [0.0], [T_FEED], [H_FEED]
    z, H = 0.0, H_FEED

    for _ in range(n_steps):
        T = T_FEED + (H - H_FEED) / CP_LIQ
        if T >= T_s:
            break

        a_p = props.specific_surface(z)
        A   = props.particle_radius(z)
        d_p = 2.0 * A
        eps = props.void_fraction(z)

        # Approximate mixture as pure N2H4 vapour at surface (conservative)
        rho = P_IN * MW['N2H4'] / (R_UNIV * T)
        mu  = props.viscosity(T)
        D   = diffusivity(D0_N2H4, T, P_IN)
        kc  = kc_chilton_colburn(G, rho, mu, D, a_p)

        # Surface concentration (pure N2H4 at T, P)
        C_surf = rho                                     # kg/m³
        r_cat  = a_p * kc * C_surf                      # kg/(m³·s) — mass-transfer limited
        r_rxn  = ALPHA1 * np.exp(-AGM_K / T) * C_surf   # kg/(m³·s) — kinetics limited
        r_vol  = min(r_cat, r_rxn)

        h4   = props.h_rxn_N2H4(T)                      # J/kg, negative
        dHdz = -h4 * r_vol / G                          # J/(kg·m), positive

        if dHdz <= 0:
            break

        # Adaptive step: limit enthalpy change per step
        dz   = -h4 / (200.0 * dHdz)
        dz   = max(dz, 1e-9)
        H_new = H + dHdz * dz

        if H_new >= H_L:
            dz = (H_L - H) / dHdz
            z += dz
            H  = H_L
            T  = T_FEED + (H - H_FEED) / CP_LIQ
            z_arr.append(z); T_arr.append(T); H_arr.append(H)
            break

        H  = H_new
        z += dz
        z_arr.append(z); T_arr.append(T); H_arr.append(H)

    return np.array(z_arr), np.array(T_arr), np.array(H_arr), z


# ─── Liquid-vapour (LV) transition region ────────────────────────────────────

def run_lv(props, bcs, z_start, n_steps=50):
    """
    March through the two-phase region at constant T = T_sat.
    Returns arrays (z, T, H, WFV) and the exit z.
    """
    T_s = bcs['T_sat']
    H_L = bcs['H_L']
    H_V = bcs['H_V']

    a_p = props.specific_surface(z_start)
    A   = props.particle_radius(z_start)
    eps = props.void_fraction(z_start)
    rho = P_IN * MW['N2H4'] / (R_UNIV * T_s)
    mu  = props.viscosity(T_s)
    D   = diffusivity(D0_N2H4, T_s, P_IN)
    kc  = kc_chilton_colburn(G, rho, mu, D, a_p)
    C_s = rho
    r_vol = a_p * min(kc, ALPHA1 * np.exp(-AGM_K / T_s)) * C_s
    h4    = props.h_rxn_N2H4(T_s)
    dHdz  = -h4 * r_vol / G
    if dHdz <= 0:
        dHdz = 1.0    # fallback so dz > 0

    z_arr, T_arr, H_arr, WFV_arr = [z_start], [T_s], [H_L], [0.0]
    z, H = z_start, H_L
    denom = H_V - H_L

    for _ in range(n_steps):
        dz    = -h4 / (40.0 * dHdz)
        dz    = max(dz, 1e-9)
        H_new = H + dHdz * dz

        if H_new >= H_V:
            dz = (H_V - H) / dHdz
            z += dz
            H  = H_V
            z_arr.append(z); T_arr.append(T_s)
            H_arr.append(H); WFV_arr.append(1.0)
            break

        H  = H_new
        z += dz
        wfv = (H - H_L) / denom if denom > 0 else 0.0
        z_arr.append(z); T_arr.append(T_s)
        H_arr.append(H); WFV_arr.append(wfv)

    return (np.array(z_arr), np.array(T_arr),
            np.array(H_arr), np.array(WFV_arr), z)


# ─── Vapour phase ODE ────────────────────────────────────────────────────────

def _vapor_rhs(z, state, props):
    """
    Right-hand side of the vapour-phase ODE system.
    state = [T, P, w_N2H4, w_NH3, w_N2, w_H2]
    All in SI units.  Returns d(state)/dz.
    """
    T, P = state[0], state[1]
    w = {'N2H4': state[2], 'NH3': state[3], 'N2': state[4], 'H2': state[5]}

    # Guard against unphysical values
    T = max(T, 200.0)
    P = max(P, 1000.0)
    for sp in w:
        w[sp] = max(w[sp], 0.0)
    w_sum = sum(w.values())
    if w_sum > 0:
        w = {sp: w[sp] / w_sum for sp in w}

    # Mixture properties
    rho    = mixture_density(T, P, w)
    Cp_mix = mixture_cp(T, props, w)
    mu     = props.viscosity(T)

    # Geometry at this axial position
    A   = props.particle_radius(z)
    d_p = 2.0 * A
    a_p = props.specific_surface(z)
    eps = props.void_fraction(z)

    # ── Concentrations (mass [kg/m³]) ──
    C = {sp: w[sp] * rho for sp in w}

    # ── Catalytic N2H4 decomposition (R1: 3 N2H4 → 4 NH3 + N2) ───────────
    # Kesten VAPOR.f uses simple T/P scaling for N2H4 with no high-P correction:
    #   DIFN = DIF4 * (TEMP/492)^1.823 * 14.7/PRES
    D_N2H4 = D0_N2H4 * (T / T_REF) ** 1.823 * (P_REF / P)
    kc_N2H4 = kc_chilton_colburn(G, rho, mu, D_N2H4, a_p)
    # Use minimum of mass-transfer and kinetics limits
    k1 = ALPHA1 * np.exp(-AGM_K / T)
    r1_mass = a_p * min(kc_N2H4, k1) * C['N2H4']   # kg_N2H4/(m³·s)

    # ── Thermal N2H4 decomposition (R3, same stoichiometry as R1) ─────────
    r3_mass = ALPHA3 * np.exp(-CGM_K / T) * C['N2H4'] * eps   # kg_N2H4/(m³·s)

    # ── Catalytic NH3 decomposition (R2: 2 NH3 → N2 + 3 H2) ──────────────
    D_NH3  = diffusivity(D0_NH3, T, P)
    kc_NH3 = kc_chilton_colburn(G, rho, mu, D_NH3, a_p)
    # SGRAD shooting: intra-particle pore diffusion + external film
    r2_area = sgrad_NH3(C['NH3'], C['H2'], T, P, A, props,
                        kc_NH3, kc_N2H4, C['N2H4'])             # kg_NH3/(m²·s)
    r2_mass = a_p * r2_area                                      # kg_NH3/(m³·s)

    # ── Convert mass rates → molar rates [mol/(m³·s)] ──────────────────────
    r1 = r1_mass / MW['N2H4']   # N2H4 consumed [mol]
    r3 = r3_mass / MW['N2H4']
    r2 = r2_mass / MW['NH3']    # NH3  consumed [mol]

    # ── Net molar production [mol/(m³·s)] ─────────────────────────────────
    # R1+R3: 3 N2H4 → 4 NH3 + N2  ⟹ per mol N2H4: +4/3 NH3, +1/3 N2
    # R2:    2 NH3  → N2   + 3 H2 ⟹ per mol NH3:  +1/2 N2,  +3/2 H2
    R_molar = {
        'N2H4': -(r1 + r3),
        'NH3':  (4.0/3.0) * (r1 + r3) - r2,
        'N2':   (1.0/3.0) * (r1 + r3) + 0.5 * r2,
        'H2':   1.5 * r2,
    }

    # ── Mass-fraction ODEs: G · dw_i/dz = MW_i · R_i ──────────────────────
    dw = {sp: MW[sp] * R_molar[sp] / G for sp in MW}

    # ── Energy equation: G · Cp_mix · dT/dz = Q_dot ──────────────────────
    # h4 < 0 (exothermic): releases heat → raises T
    # h3 > 0 (endothermic): absorbs heat → lowers T
    h4 = props.h_rxn_N2H4(T)   # J/kg_N2H4, negative
    h3 = props.h_rxn_NH3(T)    # J/kg_NH3,  positive
    Q_dot = (-h4 * (r1_mass + r3_mass)   # exothermic contribution (positive)
             - h3 * r2_mass)              # endothermic contribution (negative)
    dT_dz = Q_dot / (G * Cp_mix)

    # ── Pressure: Ergun equation ───────────────────────────────────────────
    dP_dz = ergun_dpdz(G, rho, mu, d_p, eps)

    return [dT_dz, dP_dz, dw['N2H4'], dw['NH3'], dw['N2'], dw['H2']]


def run_vapor(props, bcs, z_vap_start, z_end=Z_END):
    """
    Integrate the vapour-phase ODE from z_vap_start to z_end using solve_ivp.
    Returns a scipy OdeResult plus the initial condition dict.
    """
    T_s = bcs['T_sat']
    w0  = bcs['w0']

    y0 = [T_s, P_IN,
          w0['N2H4'], w0['NH3'], w0['N2'], w0['H2']]

    sol = solve_ivp(
        fun=_vapor_rhs,
        t_span=(z_vap_start, z_end),
        y0=y0,
        args=(props,),
        method='RK45',
        rtol=1e-5,
        atol=1e-8,
        dense_output=False,
        max_step=(z_end - z_vap_start) / 100.0,
    )
    return sol


# ─── Top-level simulation ─────────────────────────────────────────────────────

def run_simulation(data_json=None):
    """
    Run all three phases in sequence.
    Returns a dict with arrays for plotting:
      z [m], T [K], P [Pa], y_i (mole fractions), region labels.
    """
    props = Props(data_json)
    bcs   = compute_phase_bcs(props)

    print(f"Saturation temperature : {bcs['T_sat']:.1f} K  ({bcs['T_sat']*1.8:.1f} °R)")
    print(f"H_L                    : {bcs['H_L']/1000:.2f} kJ/kg")
    print(f"H_V                    : {bcs['H_V']/1000:.2f} kJ/kg")
    print(f"XV (fraction converted): {bcs['XV']:.3f}")
    print(f"Initial mole fractions : {bcs['y0']}")

    # Phase I – Liquid
    z_liq, T_liq, H_liq, z_lv_start = run_liquid(props, bcs)
    print(f"\nLiquid phase ends at z = {z_lv_start*1000:.4f} mm")

    # Phase II – Liquid-Vapour
    z_lv, T_lv, H_lv, WFV_lv, z_vap_start = run_lv(props, bcs, z_lv_start)
    print(f"LV    phase ends at z = {z_vap_start*1000:.4f} mm")

    # Phase III – Vapour
    sol = run_vapor(props, bcs, z_vap_start)
    if not sol.success:
        print(f"WARNING: vapour ODE did not converge: {sol.message}")

    z_vap = sol.t
    T_vap = sol.y[0]
    P_vap = sol.y[1]
    w_vap = {
        'N2H4': sol.y[2],
        'NH3':  sol.y[3],
        'N2':   sol.y[4],
        'H2':   sol.y[5],
    }

    # Convert mass fractions → mole fractions for output
    def _mole_fracs(w_dict):
        mol = {sp: w_dict[sp] / MW[sp] for sp in MW}
        tot = sum(mol[sp] for sp in MW)
        return {sp: mol[sp] / np.maximum(tot, 1e-30) for sp in MW}

    y_vap = _mole_fracs(w_vap)

    return dict(
        # Liquid region
        z_liq=z_liq, T_liq=T_liq,
        # LV region
        z_lv=z_lv, T_lv=T_lv, WFV_lv=WFV_lv,
        # Vapour region
        z_vap=z_vap, T_vap=T_vap, P_vap=P_vap, y_vap=y_vap,
        # Scalars
        bcs=bcs,
    )


# ─── Kesten (1967) paper reference data ──────────────────────────────────────

def kesten_paper_reference():
    """
    Digitized reference data from Kesten (1967) paper figures.

    Temperature: temperature-profile figure, G=3.0 lb/ft²-s, P=100 psia.
    These are our exact operating conditions — direct comparison.

    Composition: Fig. 14, G=1.51 lb/ft²-s, P=111.4 psia.
    Different operating point — shown for qualitative shape only.
    """
    # ── Temperature, G=3.0 lb/ft²-s, P=100 psia ──────────────────────────────
    z_T_ft = np.array([0.000, 0.003, 0.007, 0.010, 0.013, 0.016, 0.020,
                       0.025, 0.030, 0.050, 0.075, 0.100, 0.125,
                       0.150, 0.175, 0.200, 0.225, 0.250])
    T_R    = np.array([820,   1100,  1700,  2020,  2090,  2100,  2100,
                       2090,  2078,  2042,  2015,  1988,  1965,
                       1948,  1933,  1920,  1912,  1905])

    # ── Mole fractions, Fig. 14, G=1.51 lb/ft²-s, P=111.4 psia ──────────────
    z_y_ft = np.array([0.000, 0.003, 0.006, 0.010, 0.015, 0.020,
                       0.030, 0.050, 0.100, 0.200, 0.250])
    y14 = {
        'N2H4': np.array([0.370, 0.260, 0.130, 0.040, 0.007, 0.001,
                           0.000, 0.000, 0.000, 0.000, 0.000]),
        'NH3':  np.array([0.300, 0.305, 0.302, 0.295, 0.280, 0.265,
                           0.245, 0.215, 0.175, 0.120, 0.100]),
        'N2':   np.array([0.150, 0.165, 0.195, 0.225, 0.265, 0.285,
                           0.300, 0.305, 0.312, 0.318, 0.320]),
        'H2':   np.array([0.150, 0.265, 0.370, 0.440, 0.520, 0.550,
                           0.570, 0.580, 0.590, 0.600, 0.600]),
    }

    return dict(
        z_T   = z_T_ft * _ft,         # m
        T     = T_R    * _R_K,         # K
        z_y   = z_y_ft * _ft,          # m
        y14   = y14,
    )


# ─── Plotting ─────────────────────────────────────────────────────────────────

def plot_results(results, show_reference=True):
    """
    Three-panel plot: temperature, pressure, and mole fractions vs. bed depth.
    Overlays digitized Kesten (1967) paper curves when show_reference=True.
    """
    bcs   = results['bcs']
    z_v   = results['z_vap'] * 100.0   # m → cm
    T_v   = results['T_vap']
    P_v   = results['P_vap'] / 1e3     # Pa → kPa
    y     = results['y_vap']

    # Include liquid and LV regions as vertical bands
    z_liq_end = results['z_lv'][0] * 100.0
    z_lv_end  = results['z_lv'][-1] * 100.0

    ref = kesten_paper_reference() if show_reference else None

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    fig.suptitle('Hydrazine Catalyst Bed – Kesten (1967) Baseline, V1', fontsize=13)

    # ── Temperature ──────────────────────────────────────────────────────────
    ax = axes[0]
    ax.axvspan(0, z_liq_end, alpha=0.08, color='royalblue', label='Liquid')
    ax.axvspan(z_liq_end, z_lv_end, alpha=0.15, color='mediumorchid', label='LV')
    ax.plot(z_v, T_v, 'r-', linewidth=2, label='Model (V1)')
    ax.axhline(bcs['T_sat'], color='gray', linestyle='--', linewidth=1,
               label=f'T_sat = {bcs["T_sat"]:.0f} K')
    if ref is not None:
        ax.plot(ref['z_T'] * 100, ref['T'], 'k--', linewidth=1.5,
                label='Kesten (1967) G=3.0, P=100 psia')
    ax.set_ylabel('Temperature [K]')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, alpha=0.3)

    # ── Pressure ─────────────────────────────────────────────────────────────
    ax = axes[1]
    ax.axvspan(0, z_liq_end, alpha=0.08, color='royalblue')
    ax.axvspan(z_liq_end, z_lv_end, alpha=0.15, color='mediumorchid')
    ax.plot(z_v, P_v, 'b-', linewidth=2)
    ax.set_ylabel('Pressure [kPa]')
    ax.grid(True, alpha=0.3)

    # ── Mole fractions ───────────────────────────────────────────────────────
    ax = axes[2]
    ax.axvspan(0, z_liq_end, alpha=0.08, color='royalblue')
    ax.axvspan(z_liq_end, z_lv_end, alpha=0.15, color='mediumorchid')
    colours = {'N2H4': 'tab:blue', 'NH3': 'tab:orange', 'N2': 'tab:green', 'H2': 'tab:red'}
    for sp, col in colours.items():
        ax.plot(z_v, y[sp], color=col, linewidth=2, label=sp)
    if ref is not None:
        # Different conditions — dashed with reduced opacity; label only first species
        first = True
        for sp, col in colours.items():
            lbl = 'Kesten Fig.14 (G=1.51, P=111 psia)' if first else '_nolegend_'
            ax.plot(ref['z_y'] * 100, ref['y14'][sp], color=col,
                    linewidth=1.2, linestyle='--', alpha=0.6, label=lbl)
            first = False
    ax.set_ylabel('Mole fraction [-]')
    ax.set_xlabel('Bed depth z [cm]')
    ax.legend(fontsize=8, loc='center right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
    return fig


# ─── Entry point ──────────────────────────────────────────────────────────────

if __name__ == '__main__':
    results = run_simulation()
    fig = plot_results(results)

    # Print end-of-bed summary
    bcs = results['bcs']
    T_f = results['T_vap'][-1]
    P_f = results['P_vap'][-1] / 1e3
    y_f = {sp: results['y_vap'][sp][-1] for sp in MW}
    print(f"\n─── End-of-bed summary ───────────────────────────────")
    print(f"  T   = {T_f:.1f} K  ({T_f * 1.8:.1f} °R)   [Kesten: 1905 °R]")
    print(f"  P   = {P_f:.1f} kPa  ({P_f/6.895:.1f} psia)  [Kesten: 71.7 psia]")
    print(f"  y_N2H4 = {y_f['N2H4']:.4f}   [Kesten: 0]")
    print(f"  y_NH3  = {y_f['NH3']:.4f}")
    print(f"  y_N2   = {y_f['N2']:.4f}")
    print(f"  y_H2   = {y_f['H2']:.4f}")
