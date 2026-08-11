"""Reaction mechanisms for catalytic monopropellant decomposition.

Design intent
-------------
The bed solver should not know what propellant it is running. It asks a
:class:`Mechanism` for two things at every station -- species production rates
and heat release -- and everything propellant-specific lives behind that
interface. Adding hydrogen peroxide or a HAN-based propellant means writing a
new subclass, not touching the solver.

Two decisions worth stating, because they are what keep this honest:

**Heats of reaction are computed, never tabulated.** They fall out of the
species enthalpies in :mod:`fvm.chem`. A hard-coded dH can silently disagree
with the thermodynamic data it sits beside; a computed one cannot.

**Rates combine kinetics and mass transfer in series.** For a catalytic bed
the reactant must reach the surface before it can react::

    1/r = 1/r_kinetic + 1/r_mass_transfer

This matters more than it looks. Hydrazine decomposition over Shell 405 is
fast enough that it is usually diffusion-limited rather than kinetically
limited, so the answer is set by the Sherwood correlation -- which is
reasonably well established -- rather than by Arrhenius constants, which are
scattered across the literature and are the least trustworthy numbers in this
module. Where kinetics does control (ammonia dissociation), the parameters are
flagged as calibration targets rather than presented as fact.
"""
import numpy as np

from . import chem

RU = chem.RU


# ---------------------------------------------------------------------------
# rate laws
# ---------------------------------------------------------------------------
class ArrheniusRate:
    """``r = A exp(-Ea/(Ru T)) C^order``, in kmol/(m^3 bed * s).

    ``T`` is the *catalyst* temperature -- the reaction happens on the solid,
    and in a two-temperature bed the solid runs hotter than the gas in the
    decomposition zone.
    """

    def __init__(self, A, Ea, order=1.0, name=""):
        self.A = float(A)
        self.Ea = float(Ea)          # J/kmol
        self.order = float(order)
        self.name = name

    def __call__(self, T_solid, C, **_):
        C = np.maximum(C, 0.0)
        return self.A * np.exp(-self.Ea / (RU * np.maximum(T_solid, 1.0))) * C ** self.order

    def __repr__(self):
        return (f"ArrheniusRate(A={self.A:.3e}, "
                f"Ea={self.Ea / 1e6:.1f} MJ/kmol, n={self.order})")


class MassTransferRate:
    """Diffusion of reactant from the bulk gas to the catalyst surface.

    ``r = k_m a_v C``, with the mass-transfer coefficient from the
    Wakao-Funazkri packed-bed correlation

        Sh = 2 + 1.1 Sc^(1/3) Re_p^0.6

    Binary diffusivities are replaced by a fixed Schmidt number (~0.7 for
    light gases), which is accurate enough given everything else here and
    avoids carrying a diffusion database.
    """

    def __init__(self, Sc=0.7):
        self.Sc = float(Sc)

    def __call__(self, C, rho, mu, G, bed, **_):
        d_p = bed.d_p
        Re_p = np.maximum(G * d_p / np.maximum(mu, 1e-12), 1e-6)
        Sh = 2.0 + 1.1 * self.Sc ** (1.0 / 3.0) * Re_p ** 0.6
        D = mu / (np.maximum(rho, 1e-12) * self.Sc)
        k_m = Sh * D / d_p
        return k_m * bed.a_v * np.maximum(C, 0.0)


def series(r_kin, r_mt):
    """Combine kinetic and mass-transfer rates as resistances in series.

    Smooth by construction, unlike ``min(...)``, which would put a kink in the
    right-hand side and upset a stiff integrator.
    """
    r_kin = np.maximum(r_kin, 0.0)
    r_mt = np.maximum(r_mt, 0.0)
    return r_kin * r_mt / np.maximum(r_kin + r_mt, 1e-300)


# ---------------------------------------------------------------------------
# mechanism
# ---------------------------------------------------------------------------
class Reaction:
    """One reaction: stoichiometry plus a rate law.

    Stoichiometry is given as dicts of species name -> moles, written per mole
    of the *limiting reactant* so that the rate has an unambiguous meaning.
    """

    def __init__(self, reactants, products, limiting, rate_kinetic,
                 mass_transfer=None, name=""):
        self.reactants = dict(reactants)
        self.products = dict(products)
        self.limiting = limiting
        self.rate_kinetic = rate_kinetic
        self.mass_transfer = mass_transfer
        self.name = name or f"{'+'.join(self.reactants)} -> {'+'.join(self.products)}"

    def stoichiometry(self, mixture):
        """Net moles of each species per mole of reaction, as a vector."""
        nu = np.zeros(mixture.n)
        for k, v in self.reactants.items():
            nu[mixture.index(k)] -= v
        for k, v in self.products.items():
            nu[mixture.index(k)] += v
        return nu

    def check_atom_balance(self, atoms):
        """Verify element conservation. ``atoms`` maps species -> {element: n}."""
        bal = {}
        for name, v in self.reactants.items():
            for el, n in atoms[name].items():
                bal[el] = bal.get(el, 0.0) - v * n
        for name, v in self.products.items():
            for el, n in atoms[name].items():
                bal[el] = bal.get(el, 0.0) + v * n
        return {el: b for el, b in bal.items() if abs(b) > 1e-9}


class Mechanism:
    """Base class. Subclass this to add a propellant."""

    #: species participating, in a fixed order
    species_names = ()

    def __init__(self):
        self.mixture = chem.mixture(self.species_names)
        # Elemental composition comes from the same table the molecular
        # weights do, so the atom-balance check cannot drift out of step with
        # the thermodynamic data it is meant to police.
        self.atoms = {n: chem.composition(n) for n in self.species_names}
        self.reactions = self._build_reactions()
        self.nu = np.array([r.stoichiometry(self.mixture) for r in self.reactions])
        self._MW = self.mixture.MW_k

    def _build_reactions(self):
        raise NotImplementedError

    # -- inlet ------------------------------------------------------------
    def inlet_composition(self):
        """Mass fractions entering the bed."""
        raise NotImplementedError

    def inlet_enthalpy(self, T_feed):
        """Specific enthalpy of the feed [J/kg], including any phase change."""
        raise NotImplementedError

    # -- rates ------------------------------------------------------------
    def rates(self, T_gas, T_solid, p, Y, G, bed, mu=None):
        """Reaction rates [kmol/(m^3 bed * s)], one per reaction."""
        Y = np.asarray(Y, dtype=float)
        R = self.mixture.R(Y)
        rho = p / (R * np.maximum(T_gas, 1.0))
        if mu is None:
            mu = self.viscosity(T_gas, Y)
        out = []
        for rxn in self.reactions:
            i = self.mixture.index(rxn.limiting)
            C = rho * Y[i] / self._MW[i]                   # kmol/m^3
            r_k = rxn.rate_kinetic(T_solid=T_solid, C=C)
            if rxn.mass_transfer is None:
                out.append(r_k)
            else:
                r_m = rxn.mass_transfer(C=C, rho=rho, mu=mu, G=G, bed=bed)
                out.append(series(r_k, r_m))
        return np.array(out)

    def production_rates(self, r):
        """Species mass production [kg/(m^3 bed * s)] from reaction rates."""
        r = np.atleast_1d(r)
        return np.tensordot(r, self.nu, axes=(0, 0)) * self._MW

    def heat_release(self, T, r):
        """Heat released [W/m^3 bed]; positive is exothermic.

        Computed from species enthalpies, so it is automatically consistent
        with the thermodynamic data rather than a separately tabulated number.
        """
        r = np.atleast_1d(r)
        h_mole = np.array([s.h_mole(T) for s in self.mixture.species])
        dH = self.nu @ h_mole                              # J/kmol of reaction
        return -float(np.sum(r * dH)) if np.ndim(r) == 1 else -np.sum(
            r * dH[:, None], axis=0)

    def reaction_enthalpy(self, T=298.15):
        """Enthalpy of each reaction [J/kmol]; negative is exothermic."""
        h_mole = np.array([s.h_mole(T) for s in self.mixture.species])
        return self.nu @ h_mole

    # -- transport --------------------------------------------------------
    def viscosity(self, T, Y=None):
        """Mixture viscosity [Pa*s], power law about a reference point.

        Crude, but the bed answer is far more sensitive to the rate model and
        the bed geometry than to a few percent in mu.
        """
        return 3.5e-5 * (np.asarray(T, dtype=float) / 1000.0) ** 0.7

    def conductivity(self, T, Y=None):
        return self.viscosity(T, Y) * self.mixture.cp(T, Y) / 0.7

    # -- diagnostics ------------------------------------------------------
    def check_atom_balance(self):
        """Return {reaction name: imbalance} for any reaction that fails."""
        bad = {}
        for rxn in self.reactions:
            imbalance = rxn.check_atom_balance(self.atoms)
            if imbalance:
                bad[rxn.name] = imbalance
        return bad

    def __repr__(self):
        return (f"{type(self).__name__}("
                f"{len(self.reactions)} reactions, {self.mixture.n} species)")


# ---------------------------------------------------------------------------
# hydrazine
# ---------------------------------------------------------------------------
#: Enthalpy of formation of *liquid* hydrazine at 298.15 K [J/kmol].
#: The bed is fed liquid, so the 44.7 MJ/kmol vaporisation enthalpy is a real
#: energy debt the decomposition has to pay before it heats anything.
H_F_N2H4_LIQUID = 50.63e6
H_VAP_N2H4 = 44.72e6


class HydrazineShell405(Mechanism):
    """Hydrazine over an iridium-on-alumina catalyst (Shell 405 and kin).

    Two steps, written per mole of limiting reactant::

        N2H4 -> 4/3 NH3 + 1/3 N2        exothermic, very fast
        NH3  -> 1/2 N2  + 3/2 H2        endothermic, kinetically controlled

    The competition between them is the whole design problem. Step 1 sets the
    temperature; step 2 then eats some of it back while lowering the molecular
    weight. Isp peaks at partial dissociation, so a bed is sized to complete
    step 1 and control how far step 2 runs.

    Rate parameters
    ---------------
    ``A1``/``Ea1`` are deliberately set fast: step 1 over an active iridium
    catalyst is diffusion-limited in practice, so the series resistance hands
    control to the Sherwood correlation and the result is insensitive to them.

    ``A2``/``Ea2`` genuinely control the answer and are **calibration targets,
    not established constants**. Published values for ammonia dissociation over
    Shell 405 vary by orders of magnitude with catalyst age, loading and
    pretreatment. Fit them to your own bed data before trusting a dissociation
    fraction from this model.
    """

    species_names = ("N2H4", "NH3", "N2", "H2")

    def __init__(self, A1=1.0e12, Ea1=40.0e6, A2=5.0e8, Ea2=140.0e6, Sc=0.7):
        self.A1, self.Ea1 = A1, Ea1
        self.A2, self.Ea2 = A2, Ea2
        self.Sc = Sc
        super().__init__()

    def _build_reactions(self):
        return [
            Reaction({"N2H4": 1.0}, {"NH3": 4.0 / 3.0, "N2": 1.0 / 3.0},
                     limiting="N2H4",
                     rate_kinetic=ArrheniusRate(self.A1, self.Ea1, 1.0, "N2H4 decomp"),
                     mass_transfer=MassTransferRate(self.Sc),
                     name="N2H4 -> 4/3 NH3 + 1/3 N2"),
            Reaction({"NH3": 1.0}, {"N2": 0.5, "H2": 1.5},
                     limiting="NH3",
                     rate_kinetic=ArrheniusRate(self.A2, self.Ea2, 1.0, "NH3 dissoc"),
                     mass_transfer=MassTransferRate(self.Sc),
                     name="NH3 -> 1/2 N2 + 3/2 H2"),
        ]

    # -- feed -------------------------------------------------------------
    def inlet_composition(self):
        Y = np.zeros(self.mixture.n)
        Y[self.mixture.index("N2H4")] = 1.0
        return Y

    def inlet_enthalpy(self, T_feed=298.15):
        """Enthalpy of the liquid feed [J/kg].

        Referenced to the same scale as the gas-phase species, so the
        vaporisation enthalpy is carried explicitly rather than assumed away.
        """
        MW = self.mixture.MW_k[self.mixture.index("N2H4")]
        gas = chem.species("N2H4")
        h_liq_298 = H_F_N2H4_LIQUID
        cp_liq = 3080.0 * MW                     # J/(kmol*K), liquid hydrazine
        return (h_liq_298 + cp_liq * (T_feed - 298.15)) / MW

    # -- closed-form limits, useful for checking the solver ----------------
    def composition_at_dissociation(self, X):
        """Mass fractions after complete step 1 and a fraction ``X`` of step 2.

        Per mole of N2H4: ``4/3 (1-X)`` NH3, ``1/3 + 2/3 X`` N2, ``2X`` H2.
        At ``X = 1`` this is N2H4 -> N2 + 2 H2, as it must be.
        """
        X = float(X)
        return self.mixture.from_moles({
            "N2H4": 0.0,
            "NH3": 4.0 / 3.0 * (1.0 - X),
            "N2": 1.0 / 3.0 + 2.0 / 3.0 * X,
            "H2": 2.0 * X,
        })

    def adiabatic_temperature(self, X, T_feed=298.15, T_guess=1200.0):
        """Adiabatic decomposition temperature at dissociation fraction ``X``.

        The classic hydrazine design curve, and a strong check on both the
        mechanism and the thermodynamic data: it must give roughly 1650 K at
        X = 0, falling to about 880 K at X = 1.
        """
        Y = self.composition_at_dissociation(X)
        h_in = self.inlet_enthalpy(T_feed)
        return float(self.mixture.temperature_from_h(h_in, Y, T_guess=T_guess))

    def chamber_conditions(self, X, T_feed=298.15):
        """(T, MW, gamma, R, c*) at dissociation fraction ``X``.

        This is the handoff to the nozzle solver: everything it needs to build
        a frozen-composition perfect gas for the expansion.
        """
        Y = self.composition_at_dissociation(X)
        T = self.adiabatic_temperature(X, T_feed)
        MW = float(self.mixture.MW(Y))
        g = float(self.mixture.gamma(T, Y))
        R = RU / MW
        cstar = np.sqrt(g * R * T) / (g * np.sqrt((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))))
        return dict(T=T, MW=MW, gamma=g, R=R, cstar=float(cstar), Y=Y)


#: Registry so drivers can select a propellant by name.
MECHANISMS = {"hydrazine": HydrazineShell405}


def get_mechanism(name, **kw):
    try:
        return MECHANISMS[name.lower()](**kw)
    except KeyError:
        raise ValueError(f"unknown mechanism {name!r}; "
                         f"have {sorted(MECHANISMS)}") from None
