"""Thermodynamic and transport models.

Currently implements a calorically perfect (constant-cp) gas with a
frozen composition, which matches the assumptions of the algebraic
sizing tool in ``thruster_sizing.py``.
"""
import numpy as np

R_UNIVERSAL = 8314.462618  # [J/(kmol*K)]  == [J/(kg*K)] * [g/mol]
G0 = 9.80665               # [m/s^2] standard gravity


class PerfectGas:
    """Calorically perfect gas with a power-law or Sutherland viscosity.

    Parameters
    ----------
    gamma : float
        Ratio of specific heats.
    MW : float
        Molecular weight [g/mol].
    Pr : float
        Prandtl number, used to derive thermal conductivity from viscosity.
    mu_ref, T_mu_ref : float
        Reference viscosity [Pa*s] and the temperature [K] at which it applies.
    mu_law : {'power', 'sutherland'}
        ``power``      -> mu = mu_ref * (T/T_mu_ref)**omega
        ``sutherland`` -> mu = mu_ref * (T/T_mu_ref)**1.5
                               * (T_mu_ref + S_mu) / (T + S_mu)
    omega : float
        Exponent for the power law. 0.7 is representative of H2/H2O
        combustion products over a wide temperature range.
    S_mu : float
        Sutherland constant [K], only used when ``mu_law='sutherland'``.
    """

    def __init__(self, gamma=1.26, MW=11.8, Pr=0.6,
                 mu_ref=1.0e-4, T_mu_ref=3250.0,
                 mu_law="power", omega=0.7, S_mu=1064.0):
        self.gamma = float(gamma)
        self.MW = float(MW)
        self.R = R_UNIVERSAL / self.MW          # [J/(kg*K)]
        self.cv = self.R / (self.gamma - 1.0)
        self.cp = self.gamma * self.cv
        self.Pr = float(Pr)

        self.mu_ref = float(mu_ref)
        self.T_mu_ref = float(T_mu_ref)
        self.mu_law = mu_law
        self.omega = float(omega)
        self.S_mu = float(S_mu)

    # -- equation of state ------------------------------------------------
    def pressure(self, rho, e_int):
        """Static pressure from density and specific internal energy."""
        return (self.gamma - 1.0) * rho * e_int

    def temperature(self, rho, p):
        return p / (rho * self.R)

    def sound_speed(self, rho, p):
        return np.sqrt(self.gamma * p / rho)

    def enthalpy(self, rho, p, u, v):
        """Total (stagnation) enthalpy per unit mass."""
        return self.gamma / (self.gamma - 1.0) * p / rho + 0.5 * (u * u + v * v)

    # -- transport --------------------------------------------------------
    def viscosity(self, T):
        if self.mu_law == "power":
            return self.mu_ref * (T / self.T_mu_ref) ** self.omega
        elif self.mu_law == "sutherland":
            return (self.mu_ref * (T / self.T_mu_ref) ** 1.5
                    * (self.T_mu_ref + self.S_mu) / (T + self.S_mu))
        raise ValueError(f"unknown mu_law {self.mu_law!r}")

    def conductivity(self, T):
        return self.cp * self.viscosity(T) / self.Pr

    # -- isentropic relations (used for BCs, ICs and reference solutions) --
    def T_ratio(self, M):
        """T0/T."""
        return 1.0 + 0.5 * (self.gamma - 1.0) * M * M

    def p_ratio(self, M):
        """p0/p."""
        return self.T_ratio(M) ** (self.gamma / (self.gamma - 1.0))

    def rho_ratio(self, M):
        """rho0/rho."""
        return self.T_ratio(M) ** (1.0 / (self.gamma - 1.0))

    def area_ratio(self, M):
        """A/A* from the isentropic area-Mach relation."""
        g = self.gamma
        return (1.0 / M) * ((2.0 / (g + 1.0)) * self.T_ratio(M)) ** (
            (g + 1.0) / (2.0 * (g - 1.0)))

    def cstar_ideal(self, T0):
        """Ideal characteristic velocity [m/s]."""
        g = self.gamma
        return np.sqrt(g * self.R * T0) / (
            g * np.sqrt((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))))

    def __repr__(self):
        return (f"PerfectGas(gamma={self.gamma}, MW={self.MW}, "
                f"R={self.R:.2f} J/kg/K, cp={self.cp:.1f} J/kg/K, Pr={self.Pr})")


#: LOX/LH2 at O/F = 5.0 -- the propellant combination used by thruster_sizing.py
LOX_LH2_OF5 = dict(gamma=1.26, MW=11.8, Pr=0.6,
                   mu_ref=1.0e-4, T_mu_ref=3250.0, mu_law="power", omega=0.7)

#: Cold gaseous nitrogen. This is the working fluid in Hayn's micronozzle
#: experiments (NASA TM-77730), so it is what the validation case needs.
#: Sutherland rather than a power law, because the temperature range is
#: modest and Sutherland is accurate for N2 near ambient.
COLD_N2 = dict(gamma=1.4, MW=28.0134, Pr=0.72,
               mu_ref=1.78e-5, T_mu_ref=300.0, mu_law="sutherland", S_mu=111.0)
