"""Reaction mechanisms for catalytic monopropellant decomposition.

Design intent
-------------
The bed solver should not know what propellant it is running. It asks a
:class:`Mechanism` for species production rates and heat release, and
everything propellant-specific lives behind that interface. Adding hydrogen
peroxide or a HAN-based propellant means writing a new subclass, not touching
the solver.

Heats of reaction are computed from the species enthalpies in :mod:`fvm.chem`,
never tabulated. A hard-coded dH can silently disagree with the thermodynamic
data it sits beside; a derived one cannot.

Provenance of the hydrazine model
---------------------------------
Rate parameters and correlations come from Kesten's UARL work under NASA
contract NAS 7-458, both reports being in ``docs/``:

* F910461-12, *Analytical Study of Catalytic Reactors for Hydrazine
  Decomposition*, First Annual Progress Report, May 1967.
* G910461-30, *Computer Programs Manual*, August 1968.

Kesten is candid about what these numbers are, and so is this module:

* The hydrazine catalytic activation energy (2500 degR) was, in his words,
  "chosen rather arbitrarily"; the pre-exponential was then fitted to engine
  test data. This matters less than it sounds, because that step is
  diffusion-controlled -- a claim his own report makes and which a test here
  asserts.
* The ammonia rate's hydrogen-inhibition order was measured on *platinum*
  (Melton; Logan and Kemball) and assumed to carry over to Shell 405.
  Kesten: "this assumption remains untested".
* The ammonia pre-exponential is a fit, quoted as 0.3e11 from his steady-state
  model and 1e11 from his transient model, with the true value "probably
  between" them. That factor of three is the honest uncertainty on the single
  parameter that most controls predicted dissociation.

Treat all of them as well-founded defaults to be re-fitted against your own
bed data, not as physical constants.
"""
import numpy as np

from . import chem

RU = chem.RU

# -- unit conversions to and from Kesten's units ----------------------------
LB_FT3_TO_KG_M3 = 16.018463374
FT2_S_TO_M2_S = 0.09290304
RANKINE_TO_KELVIN = 5.0 / 9.0
PSIA_TO_PA = 6894.757293


def activation_from_degR(Ea_over_R_degR):
    """Convert Kesten's ``Ea/R`` in deg R to an activation energy in J/kmol."""
    return RU * Ea_over_R_degR * RANKINE_TO_KELVIN


def prefactor_to_SI(A_imperial, net_concentration_order):
    """Convert a rate pre-exponential from Kesten's units into SI.

    First-order rates need no conversion: ``r = A C`` has ``A`` in 1/s in any
    unit system. Rates whose concentration exponents do not sum to one do
    need it -- the ammonia rate is first order in NH3 and order -1.6 in H2, so
    its pre-exponential carries units of ``(lb/ft^3)^1.6 / s`` and converting
    it is not optional.

    ``net_concentration_order`` is the sum of all concentration exponents
    (1.0 - 1.6 = -0.6 for the ammonia rate).
    """
    return A_imperial * LB_FT3_TO_KG_M3 ** (1.0 - net_concentration_order)


class RateState:
    """Everything a rate law may need at one station."""

    __slots__ = ("T_gas", "T_solid", "p", "C", "rho", "mu", "G", "bed", "mixture")

    def __init__(self, T_gas, T_solid, p, C, rho, mu, G, bed, mixture):
        self.T_gas, self.T_solid, self.p = T_gas, T_solid, p
        self.C, self.rho, self.mu = C, rho, mu
        self.G, self.bed, self.mixture = G, bed, mixture

    def conc(self, name):
        """Mass concentration of a species [kg/m^3]."""
        return self.C[self.mixture.index(name)]


# ---------------------------------------------------------------------------
# rate laws
# ---------------------------------------------------------------------------
class CatalyticRate:
    """``A exp(-Ea/(Ru T_s)) C_react^m / C_inhib^n``, in kg/(m^3 bed * s).

    The temperature is the *catalyst* temperature: the reaction happens on the
    solid, which in a two-temperature bed runs hotter than the gas wherever
    decomposition is active.

    ``inhibitor`` implements product inhibition -- for ammonia over platinum
    and (by assumption) Shell 405, hydrogen suppresses the reaction, which is
    a self-limiting feedback that keeps a bed from running away toward complete
    dissociation.
    """

    def __init__(self, A, Ea, reactant, order=1.0,
                 inhibitor=None, inhibitor_order=0.0, inhibitor_floor=0.0,
                 name=""):
        self.A = float(A)
        self.Ea = float(Ea)
        self.reactant = reactant
        self.order = float(order)
        self.inhibitor = inhibitor
        self.inhibitor_order = float(inhibitor_order)
        self.inhibitor_floor = float(inhibitor_floor)
        self.name = name

    def __call__(self, st):
        C = np.maximum(st.conc(self.reactant), 0.0)
        r = self.A * np.exp(-self.Ea / (RU * np.maximum(st.T_solid, 1.0))) * C ** self.order
        if self.inhibitor is not None:
            Ci = np.maximum(st.conc(self.inhibitor), self.inhibitor_floor)
            r = r / Ci ** self.inhibitor_order
        return r

    def __repr__(self):
        s = f"CatalyticRate(A={self.A:.3e}, Ea={self.Ea / 1e6:.2f} MJ/kmol"
        if self.inhibitor:
            s += f", 1/[{self.inhibitor}]^{self.inhibitor_order}"
        return s + ")"


class HomogeneousRate(CatalyticRate):
    """Same form, but evaluated at the *gas* temperature and with no catalyst.

    Kesten carries a gas-phase thermal decomposition path alongside the
    catalytic one; it is negligible cold and contributes once the bed is hot.
    """

    def __call__(self, st):
        C = np.maximum(st.conc(self.reactant), 0.0)
        return self.A * np.exp(-self.Ea / (RU * np.maximum(st.T_gas, 1.0))) * C ** self.order


class KestenMassTransfer:
    """Film diffusion to the catalyst surface, using Kesten's own correlation.

    From subroutine ``KCF`` in the UARL listings::

        k_c = 0.61 (G/rho) Sc^-0.667 [G/(a_p mu)]^-0.41

    a Colburn j-factor form with the Reynolds number built on surface area per
    unit bed volume rather than particle diameter. Every group is
    dimensionless except ``G/rho``, so the correlation carries across unit
    systems unchanged.

    Diffusivity is corrected from its STP value the way the Fortran does::

        D(T, p) = D_STP (p_ref/p) (T/T_ref)^1.823
    """

    T_REF = 492.0 * RANKINE_TO_KELVIN        # 273.15 K
    P_REF = 14.7 * PSIA_TO_PA                # 101325 Pa

    def __init__(self, D_stp, species):
        self.D_stp = float(D_stp)            # m^2/s at T_REF, P_REF
        self.species = species

    def diffusivity(self, T, p):
        return self.D_stp * (self.P_REF / np.maximum(p, 1.0)) * (T / self.T_REF) ** 1.823

    def __call__(self, st):
        D = self.diffusivity(st.T_gas, st.p)
        rho = np.maximum(st.rho, 1e-12)
        mu = np.maximum(st.mu, 1e-12)
        Sc = mu / (rho * D)
        Re_a = np.maximum(st.G / (st.bed.a_v * mu), 1e-12)
        kc = 0.61 * (st.G / rho) * Sc ** -0.667 * Re_a ** -0.41
        return kc * st.bed.a_v * np.maximum(st.conc(self.species), 0.0)


def series(r_kin, r_mt):
    """Combine kinetic and mass-transfer rates as resistances in series.

    Smooth by construction, unlike ``min(...)``, which would put a kink in the
    right-hand side and upset a stiff integrator. It also bounds the rate: no
    reaction can outrun the supply of reactant to the surface, which is what
    keeps the hydrogen-inhibited ammonia rate finite as H2 gets small.
    """
    r_kin = np.maximum(r_kin, 0.0)
    r_mt = np.maximum(r_mt, 0.0)
    return r_kin * r_mt / np.maximum(r_kin + r_mt, 1e-300)


# ---------------------------------------------------------------------------
# mechanism
# ---------------------------------------------------------------------------
class Reaction:
    """One reaction: stoichiometry plus a rate law.

    Stoichiometry is written per mole of the *limiting reactant*, so the rate
    has an unambiguous meaning. Rate laws return a mass rate for that species,
    which is divided by its molecular weight to give reaction extent.
    """

    def __init__(self, reactants, products, limiting, rate,
                 mass_transfer=None, name=""):
        self.reactants = dict(reactants)
        self.products = dict(products)
        self.limiting = limiting
        self.rate = rate
        self.mass_transfer = mass_transfer
        self.name = name or f"{'+'.join(self.reactants)} -> {'+'.join(self.products)}"

    def stoichiometry(self, mixture):
        nu = np.zeros(mixture.n)
        for k, v in self.reactants.items():
            nu[mixture.index(k)] -= v
        for k, v in self.products.items():
            nu[mixture.index(k)] += v
        return nu

    def check_atom_balance(self, atoms):
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

    def inlet_composition(self):
        raise NotImplementedError

    def inlet_enthalpy(self, T_feed):
        raise NotImplementedError

    # -- rates ------------------------------------------------------------
    def rates(self, T_gas, T_solid, p, Y, G, bed, mu=None):
        """Reaction extents [kmol/(m^3 bed * s)], one per reaction."""
        Y = np.asarray(Y, dtype=float)
        rho = p / (self.mixture.R(Y) * np.maximum(T_gas, 1.0))
        if mu is None:
            mu = self.viscosity(T_gas, Y)
        st = RateState(T_gas, T_solid, p, rho * Y, rho, mu, G, bed, self.mixture)

        out = []
        for rxn in self.reactions:
            r = rxn.rate(st)                                   # kg/(m^3 s)
            if rxn.mass_transfer is not None:
                r = series(r, rxn.mass_transfer(st))
            out.append(r / self._MW[self.mixture.index(rxn.limiting)])
        return np.array(out)

    def production_rates(self, r):
        """Species mass production [kg/(m^3 bed * s)] from reaction extents."""
        r = np.atleast_1d(r)
        return np.tensordot(r, self.nu, axes=(0, 0)) * self._MW

    def heat_release(self, T, r):
        """Heat released [W/m^3 bed]; positive is exothermic."""
        r = np.atleast_1d(r)
        h_mole = np.array([s.h_mole(T) for s in self.mixture.species])
        dH = self.nu @ h_mole
        return -float(np.sum(r * dH)) if np.ndim(r) == 1 else -np.sum(
            r * dH[:, None], axis=0)

    def reaction_enthalpy(self, T=298.15):
        """Enthalpy of each reaction [J/kmol]; negative is exothermic."""
        h_mole = np.array([s.h_mole(T) for s in self.mixture.species])
        return self.nu @ h_mole

    # -- transport --------------------------------------------------------
    def viscosity(self, T, Y=None):
        return 3.5e-5 * (np.asarray(T, dtype=float) / 1000.0) ** 0.7

    def conductivity(self, T, Y=None):
        return self.viscosity(T, Y) * self.mixture.cp(T, Y) / 0.7

    def check_atom_balance(self):
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
#: The bed is fed liquid, so vaporisation is a real energy debt the
#: decomposition has to pay before it heats anything.
H_F_N2H4_LIQUID = 50.63e6

#: Vaporisation enthalpy is *derived*, not quoted, so that
#: ``h_f(liquid) + h_vap`` always equals the gas-phase species enthalpy the
#: reactions are written against. Quoting both independently left a 0.17
#: kJ/mol inconsistency against the CEA value of 95.18 -- small, but there is
#: no reason to carry a redundant constant that can drift.
H_VAP_N2H4 = chem.species("N2H4").h_formation() - H_F_N2H4_LIQUID

CP_N2H4_LIQUID = 0.7332 * 4186.8       # J/(kg*K) -- Kesten's CFL, 0.7332 Btu/lb-degR

#: Kesten's parameters, in his units. Converted on construction.
KESTEN = dict(
    A_N2H4_cat=1.0e10,          # 1/s
    EaR_N2H4_cat=2500.0,        # deg R   ("chosen rather arbitrarily")
    A_NH3_cat=0.3e11,           # (lb/ft^3)^1.6 / s   (steady-state fit; 1e11 transient)
    EaR_NH3_cat=50000.0,        # deg R
    A_N2H4_hom=2.14e10,         # 1/s
    EaR_N2H4_hom=33000.0,       # deg R
    D_N2H4_stp=0.95e-4,         # ft^2/s at STP
    D_NH3_stp=0.17e-3,          # ft^2/s at STP
)


class HydrazineShell405(Mechanism):
    """Hydrazine over an iridium-on-alumina catalyst (Shell 405 and kin).

    Three paths, each written per mole of its limiting reactant::

        N2H4 -> 4/3 NH3 + 1/3 N2      catalytic, diffusion-controlled
        N2H4 -> 4/3 NH3 + 1/3 N2      homogeneous, gas phase
        NH3  -> 1/2 N2  + 3/2 H2      catalytic, kinetically controlled,
                                      inhibited by hydrogen

    The competition is the design problem. Decomposition sets the temperature;
    dissociation then eats some of it back while lowering the molecular
    weight, so c* peaks at partial dissociation and a bed is sized to complete
    the first step and control how far the second runs.

    Parameters
    ----------
    nH2 : float
        Hydrogen inhibition order. Kesten's own sources disagree -- Melton
        reports 1.0, Logan and Kemball 1.6 -- and so do two copies of his
        Fortran: ``PARAM.f`` divides by ``C1**1.6`` while ``MAIN.f`` uses
        ``C1**1.0``. The accompanying ``kinetics.json`` picks 1.6, which is
        the default here.
    A_NH3 : float
        Ammonia pre-exponential in Kesten's units. Default is his
        steady-state fit; his transient model gives 1e11 and he expects the
        truth in between. This is the parameter worth fitting first.
    Y_H2_floor : float
        Hydrogen mass fraction floor used only to keep ``1/C_H2^n`` finite
        where the bed has produced ammonia but almost no hydrogen. The
        original Fortran has no such floor, so the singularity is real and
        merely never reached by its integration path. Results should be
        checked for sensitivity to this value; the series resistance with mass
        transfer provides the physical bound, this only removes the numerical
        one.
    """

    species_names = ("N2H4", "NH3", "N2", "H2")

    def __init__(self, nH2=1.6, A_NH3=None, Y_H2_floor=1e-4,
                 include_homogeneous=True, **overrides):
        k = dict(KESTEN)
        k.update(overrides)
        self.k = k
        self.nH2 = float(nH2)
        self.A_NH3_imperial = float(k["A_NH3_cat"] if A_NH3 is None else A_NH3)
        self.Y_H2_floor = float(Y_H2_floor)
        self.include_homogeneous = bool(include_homogeneous)
        super().__init__()

    def _build_reactions(self):
        k = self.k
        decomp_products = {"NH3": 4.0 / 3.0, "N2": 1.0 / 3.0}

        # The ammonia pre-exponential is the only one whose units depend on
        # the unit system: its concentration exponents sum to 1 - nH2.
        A_NH3_SI = prefactor_to_SI(self.A_NH3_imperial, 1.0 - self.nH2)
        rho_ref = 1.0                       # floor is applied as a concentration
        rxns = [
            Reaction({"N2H4": 1.0}, decomp_products, limiting="N2H4",
                     rate=CatalyticRate(k["A_N2H4_cat"],
                                        activation_from_degR(k["EaR_N2H4_cat"]),
                                        reactant="N2H4", name="N2H4 catalytic"),
                     mass_transfer=KestenMassTransfer(
                         k["D_N2H4_stp"] * FT2_S_TO_M2_S, "N2H4"),
                     name="N2H4 -> 4/3 NH3 + 1/3 N2 (catalytic)"),
            Reaction({"NH3": 1.0}, {"N2": 0.5, "H2": 1.5}, limiting="NH3",
                     rate=CatalyticRate(A_NH3_SI,
                                        activation_from_degR(k["EaR_NH3_cat"]),
                                        reactant="NH3",
                                        inhibitor="H2", inhibitor_order=self.nH2,
                                        inhibitor_floor=self.Y_H2_floor * rho_ref,
                                        name="NH3 catalytic"),
                     mass_transfer=KestenMassTransfer(
                         k["D_NH3_stp"] * FT2_S_TO_M2_S, "NH3"),
                     name="NH3 -> 1/2 N2 + 3/2 H2 (catalytic)"),
        ]
        if self.include_homogeneous:
            rxns.append(
                Reaction({"N2H4": 1.0}, decomp_products, limiting="N2H4",
                         rate=HomogeneousRate(k["A_N2H4_hom"],
                                              activation_from_degR(k["EaR_N2H4_hom"]),
                                              reactant="N2H4", name="N2H4 thermal"),
                         mass_transfer=None,
                         name="N2H4 -> 4/3 NH3 + 1/3 N2 (homogeneous)"))
        return rxns

    # -- feed -------------------------------------------------------------
    def inlet_composition(self):
        Y = np.zeros(self.mixture.n)
        Y[self.mixture.index("N2H4")] = 1.0
        return Y

    def inlet_enthalpy(self, T_feed=298.15):
        """Enthalpy of the liquid feed [J/kg], on the gas-phase scale."""
        MW = self.mixture.MW_k[self.mixture.index("N2H4")]
        return (H_F_N2H4_LIQUID + CP_N2H4_LIQUID * MW * (T_feed - 298.15)) / MW

    # -- closed-form limits, useful for checking the solver ----------------
    @staticmethod
    def X_from_kesten_f(f):
        """Convert Kesten's dissociation fraction to the two-step convention.

        His footnote (F910461-12 p.11 / G910461-30 p.5) gives
        ``f_two_step = (3 f + 1) / 4``. His experimentally determined overall
        reaction already implies 25% of the ammonia has dissociated, so his
        ``f = 0`` is ``X = 0.25``. Comparing the two directly is wrong by
        roughly 90 K at typical bed conditions.
        """
        return (3.0 * np.asarray(f, dtype=float) + 1.0) / 4.0

    @staticmethod
    def kesten_f_from_X(X):
        return (4.0 * np.asarray(X, dtype=float) - 1.0) / 3.0

    def composition_at_dissociation(self, X):
        """Mass fractions after complete step 1 and a fraction ``X`` of step 2."""
        X = float(X)
        return self.mixture.from_moles({
            "N2H4": 0.0,
            "NH3": 4.0 / 3.0 * (1.0 - X),
            "N2": 1.0 / 3.0 + 2.0 / 3.0 * X,
            "H2": 2.0 * X,
        })

    def adiabatic_temperature(self, X, T_feed=298.15, T_guess=1200.0):
        """Adiabatic decomposition temperature at dissociation fraction ``X``."""
        Y = self.composition_at_dissociation(X)
        h_in = self.inlet_enthalpy(T_feed)
        return float(self.mixture.temperature_from_h(h_in, Y, T_guess=T_guess))

    def chamber_conditions(self, X, T_feed=298.15):
        """(T, MW, gamma, R, c*) at dissociation fraction ``X``."""
        Y = self.composition_at_dissociation(X)
        T = self.adiabatic_temperature(X, T_feed)
        MW = float(self.mixture.MW(Y))
        g = float(self.mixture.gamma(T, Y))
        R = RU / MW
        cstar = np.sqrt(g * R * T) / (g * np.sqrt((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))))
        return dict(T=T, MW=MW, gamma=g, R=R, cstar=float(cstar), Y=Y)


MECHANISMS = {"hydrazine": HydrazineShell405}


def get_mechanism(name, **kw):
    try:
        return MECHANISMS[name.lower()](**kw)
    except KeyError:
        raise ValueError(f"unknown mechanism {name!r}; "
                         f"have {sorted(MECHANISMS)}") from None
