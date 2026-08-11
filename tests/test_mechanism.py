"""Tests for reaction mechanisms.

The load-bearing test here is the hydrazine adiabatic decomposition curve.
It is a well-known design curve -- roughly 1650 K with no ammonia
dissociation, falling to about 880 K when dissociation is complete -- and
reproducing it exercises the stoichiometry, the liquid-feed enthalpy including
vaporisation, and the NASA-7 thermodynamic data all at once. If any of those
is wrong the curve misses.
"""
import numpy as np
import pytest

from fvm import chem
from fvm.mechanism import (ArrheniusRate, HydrazineShell405, MassTransferRate,
                           get_mechanism, series)


class FakeBed:
    """Minimal stand-in for the bed geometry of stage C."""
    d_p = 0.3e-3        # m
    a_v = 6000.0        # m^2 per m^3 of bed
    eps = 0.4


@pytest.fixture(scope="module")
def mech():
    return HydrazineShell405()


# -- stoichiometry ----------------------------------------------------------
def test_reactions_balance_atoms(mech):
    assert mech.check_atom_balance() == {}


def test_reactions_conserve_mass(mech):
    """sum(nu_i * MW_i) must vanish for every reaction."""
    for j, rxn in enumerate(mech.reactions):
        residual = float(mech.nu[j] @ mech.mixture.MW_k)
        assert abs(residual) < 1e-9, f"{rxn.name} loses {residual:.3e} kg/kmol"


def test_production_rates_conserve_mass(mech):
    r = np.array([2.0, 0.7])
    assert abs(mech.production_rates(r).sum()) < 1e-9


def test_step_one_is_exothermic_and_step_two_endothermic(mech):
    dH = mech.reaction_enthalpy(298.15)
    assert dH[0] < 0.0, "N2H4 decomposition must release heat"
    assert dH[1] > 0.0, "NH3 dissociation must absorb heat"


def test_complete_dissociation_gives_nitrogen_and_hydrogen(mech):
    """At X = 1 the overall reaction is N2H4 -> N2 + 2 H2."""
    Y = mech.composition_at_dissociation(1.0)
    x = mech.mixture.mole_fractions(Y)
    assert x[mech.mixture.index("NH3")] == pytest.approx(0.0, abs=1e-12)
    assert x[mech.mixture.index("N2")] == pytest.approx(1.0 / 3.0, rel=1e-9)
    assert x[mech.mixture.index("H2")] == pytest.approx(2.0 / 3.0, rel=1e-9)


def test_nitrogen_to_hydrogen_ratio_is_fixed_by_the_feed(mech):
    """Every X must preserve the 2:4 N:H ratio of N2H4."""
    for X in (0.0, 0.25, 0.6, 1.0):
        Y = mech.composition_at_dissociation(X)
        x = mech.mixture.mole_fractions(Y)
        N = 2 * x[mech.mixture.index("N2")] + x[mech.mixture.index("NH3")]
        H = 2 * x[mech.mixture.index("H2")] + 3 * x[mech.mixture.index("NH3")]
        assert H / N == pytest.approx(2.0, rel=1e-9)


# -- the design curve -------------------------------------------------------
def test_adiabatic_temperature_matches_known_hydrazine_values(mech):
    """Roughly 1650 K undissociated, about 880 K fully dissociated.

    Tolerances are wide because published values themselves scatter by a few
    tens of K depending on the reference state used.
    """
    assert mech.adiabatic_temperature(0.0) == pytest.approx(1650.0, abs=90.0)
    assert mech.adiabatic_temperature(1.0) == pytest.approx(880.0, abs=70.0)


def test_adiabatic_temperature_falls_monotonically_with_dissociation(mech):
    X = np.linspace(0.0, 1.0, 11)
    T = np.array([mech.adiabatic_temperature(x) for x in X])
    assert np.all(np.diff(T) < 0.0), "NH3 dissociation is endothermic"


def test_molecular_weight_falls_with_dissociation(mech):
    X = np.linspace(0.0, 1.0, 11)
    MW = np.array([mech.chamber_conditions(x)["MW"] for x in X])
    assert np.all(np.diff(MW) < 0.0)
    assert MW[0] == pytest.approx(19.2, abs=0.5)     # 4/3 NH3 + 1/3 N2
    assert MW[-1] == pytest.approx(10.7, abs=0.3)    # N2 + 2 H2


def test_cstar_is_plausible_and_has_an_interior_optimum(mech):
    """c* trades falling temperature against falling molecular weight.

    Both effects are strong and they oppose, which is why hydrazine engines
    are designed for partial dissociation rather than either extreme.
    """
    X = np.linspace(0.0, 1.0, 21)
    cstar = np.array([mech.chamber_conditions(x)["cstar"] for x in X])
    assert np.all((cstar > 1100.0) & (cstar < 1600.0))
    assert 0 < int(np.argmax(cstar)) < len(X) - 1, "optimum should be interior"


def test_liquid_feed_costs_the_vaporisation_enthalpy(mech):
    """Feeding liquid must give a lower flame temperature than feeding gas."""
    from fvm.mechanism import H_F_N2H4_LIQUID, H_VAP_N2H4
    MW = mech.mixture.MW_k[mech.mixture.index("N2H4")]
    Y = mech.composition_at_dissociation(0.0)
    T_liq = mech.adiabatic_temperature(0.0)
    T_gas = float(mech.mixture.temperature_from_h(
        (H_F_N2H4_LIQUID + H_VAP_N2H4) / MW, Y, T_guess=1500.0))
    assert T_gas > T_liq + 100.0


# -- rate laws --------------------------------------------------------------
def test_arrhenius_increases_with_temperature_and_concentration():
    rate = ArrheniusRate(1e10, 100e6, order=1.0)
    assert rate(T_solid=1200.0, C=1.0) > rate(T_solid=800.0, C=1.0)
    assert rate(T_solid=1000.0, C=2.0) == pytest.approx(2 * rate(T_solid=1000.0, C=1.0))


def test_series_resistance_is_bounded_by_both_branches():
    a = np.array([1.0, 100.0, 1e6, 0.0])
    b = np.array([100.0, 1.0, 1.0, 5.0])
    s = series(a, b)
    assert np.all(s <= a + 1e-12)
    assert np.all(s <= b + 1e-12)


def test_series_resistance_tends_to_the_slower_branch():
    """A very fast kinetic rate must hand control to mass transfer."""
    assert series(1e12, 3.0) == pytest.approx(3.0, rel=1e-6)
    assert series(3.0, 1e12) == pytest.approx(3.0, rel=1e-6)


def test_mass_transfer_rate_rises_with_mass_flux():
    mt = MassTransferRate()
    bed = FakeBed()
    lo = mt(C=1.0, rho=2.0, mu=3e-5, G=5.0, bed=bed)
    hi = mt(C=1.0, rho=2.0, mu=3e-5, G=50.0, bed=bed)
    assert hi > lo


def test_rates_are_nonnegative_and_finite(mech):
    Y = mech.mixture.from_moles({"N2H4": 1.0})
    r = mech.rates(T_gas=600.0, T_solid=900.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    assert np.all(np.isfinite(r)) and np.all(r >= 0.0)


def test_no_reaction_without_reactant(mech):
    """Rates must vanish when the limiting species is absent."""
    Y = mech.composition_at_dissociation(1.0)     # no N2H4, no NH3 left
    r = mech.rates(T_gas=1000.0, T_solid=1000.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    assert np.allclose(r, 0.0, atol=1e-12)


def test_heat_release_is_positive_while_decomposing(mech):
    Y = mech.mixture.from_moles({"N2H4": 1.0})
    r = mech.rates(T_gas=700.0, T_solid=1000.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    assert mech.heat_release(1000.0, r) > 0.0


def test_fast_step_one_is_diffusion_limited(mech):
    """The design claim: step 1's answer must not depend on its Arrhenius A.

    If it does, the rate constants -- the least trustworthy numbers in the
    module -- are controlling the result.
    """
    Y = mech.mixture.from_moles({"N2H4": 1.0})
    kw = dict(T_gas=800.0, T_solid=1100.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    slow = HydrazineShell405(A1=1.0e11).rates(**kw)[0]
    fast = HydrazineShell405(A1=1.0e14).rates(**kw)[0]
    assert fast == pytest.approx(slow, rel=0.02)


# -- registry ---------------------------------------------------------------
def test_registry_returns_a_mechanism():
    assert isinstance(get_mechanism("hydrazine"), HydrazineShell405)


def test_registry_rejects_unknown_names():
    with pytest.raises(ValueError, match="unknown mechanism"):
        get_mechanism("unobtainium")
