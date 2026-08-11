"""Tests for reaction mechanisms.

Two load-bearing checks here:

* The hydrazine adiabatic decomposition curve -- roughly 1650 K undissociated
  falling to 880 K fully dissociated -- which exercises stoichiometry, the
  liquid-feed enthalpy including vaporisation, and the NASA-7 data together.
* Kesten's own reactor output, in ``docs/kesten_claude``, which is an
  independent implementation of the same physics from 1967. Agreeing with it
  is worth more than any amount of internal consistency.
"""
import numpy as np
import pytest

from fvm import chem
from fvm.mechanism import (CatalyticRate, HomogeneousRate, HydrazineShell405,
                           KestenMassTransfer, RateState, activation_from_degR,
                           get_mechanism, prefactor_to_SI, series,
                           LB_FT3_TO_KG_M3, RANKINE_TO_KELVIN)


class FakeBed:
    """Minimal stand-in for the bed geometry of stage C."""
    d_p = 0.3e-3        # m
    a_v = 6000.0        # m^2 per m^3 of bed
    eps = 0.4


@pytest.fixture(scope="module")
def mech():
    return HydrazineShell405()


def _state(mech, Y, T_gas=800.0, T_solid=1000.0, p=8e5, G=20.0):
    rho = p / (mech.mixture.R(Y) * T_gas)
    return RateState(T_gas, T_solid, p, rho * np.asarray(Y), rho,
                     mech.viscosity(T_gas, Y), G, FakeBed(), mech.mixture)


# -- unit conversions -------------------------------------------------------
def test_activation_energy_conversion_from_degR():
    """2500 degR is 1388.9 K, so Ea = Ru * 1388.9."""
    assert activation_from_degR(2500.0) == pytest.approx(chem.RU * 2500 * 5 / 9)
    assert activation_from_degR(50000.0) / 1e6 == pytest.approx(230.96, rel=1e-3)


def test_first_order_prefactor_needs_no_conversion():
    """r = A C has A in 1/s in any unit system."""
    assert prefactor_to_SI(1.0e10, 1.0) == pytest.approx(1.0e10)


def test_ammonia_prefactor_conversion_matches_hand_calculation():
    """Exponents sum to 1 - 1.6 = -0.6, so A scales by (lb/ft^3 -> kg/m^3)^1.6."""
    A_si = prefactor_to_SI(0.3e11, 1.0 - 1.6)
    assert A_si == pytest.approx(0.3e11 * LB_FT3_TO_KG_M3 ** 1.6, rel=1e-12)


def test_rate_is_unit_system_independent():
    """The ammonia rate must give the same physical answer in either system.

    This is the check that catches a botched pre-exponential conversion: build
    the rate in SI, evaluate it at concentrations expressed in kg/m^3, and
    compare against Kesten's expression evaluated in lb/ft^3.
    """
    n, A_imp, EaR = 1.6, 0.3e11, 50000.0
    C_NH3_lb, C_H2_lb, T_degR = 0.05, 0.004, 1900.0

    r_imperial = A_imp * np.exp(-EaR / T_degR) * C_NH3_lb / C_H2_lb ** n   # lb/ft^3-s
    r_expected_si = r_imperial * LB_FT3_TO_KG_M3                           # kg/m^3-s

    A_si = prefactor_to_SI(A_imp, 1.0 - n)
    Ea = activation_from_degR(EaR)
    r_si = (A_si * np.exp(-Ea / (chem.RU * T_degR * RANKINE_TO_KELVIN))
            * (C_NH3_lb * LB_FT3_TO_KG_M3) / (C_H2_lb * LB_FT3_TO_KG_M3) ** n)

    assert r_si == pytest.approx(r_expected_si, rel=1e-10)


# -- stoichiometry ----------------------------------------------------------
def test_reactions_balance_atoms(mech):
    assert mech.check_atom_balance() == {}


def test_reactions_conserve_mass(mech):
    for j, rxn in enumerate(mech.reactions):
        residual = float(mech.nu[j] @ mech.mixture.MW_k)
        assert abs(residual) < 1e-9, f"{rxn.name} loses {residual:.3e} kg/kmol"


def test_production_rates_conserve_mass(mech):
    r = np.array([2.0, 0.7, 0.1])
    assert abs(mech.production_rates(r).sum()) < 1e-9


def test_step_one_is_exothermic_and_step_two_endothermic(mech):
    dH = mech.reaction_enthalpy(298.15)
    assert dH[0] < 0.0
    assert dH[1] > 0.0
    assert dH[2] < 0.0, "the homogeneous path is the same reaction, so also exothermic"


def test_complete_dissociation_gives_nitrogen_and_hydrogen(mech):
    Y = mech.composition_at_dissociation(1.0)
    x = mech.mixture.mole_fractions(Y)
    assert x[mech.mixture.index("NH3")] == pytest.approx(0.0, abs=1e-12)
    assert x[mech.mixture.index("N2")] == pytest.approx(1.0 / 3.0, rel=1e-9)
    assert x[mech.mixture.index("H2")] == pytest.approx(2.0 / 3.0, rel=1e-9)


# -- the design curve -------------------------------------------------------
def test_adiabatic_temperature_matches_known_hydrazine_values(mech):
    assert mech.adiabatic_temperature(0.0) == pytest.approx(1650.0, abs=90.0)
    assert mech.adiabatic_temperature(1.0) == pytest.approx(880.0, abs=70.0)


def test_adiabatic_temperature_falls_monotonically_with_dissociation(mech):
    T = np.array([mech.adiabatic_temperature(x) for x in np.linspace(0, 1, 11)])
    assert np.all(np.diff(T) < 0.0)


def test_cstar_has_an_interior_optimum(mech):
    X = np.linspace(0.0, 1.0, 21)
    cstar = np.array([mech.chamber_conditions(x)["cstar"] for x in X])
    assert np.all((cstar > 1100.0) & (cstar < 1600.0))
    assert 0 < int(np.argmax(cstar)) < len(X) - 1


# -- Kesten's dissociation-fraction convention ------------------------------
def test_kesten_f_of_zero_is_a_quarter_dissociated(mech):
    """His overall reaction already implies 25% of the ammonia has gone."""
    assert mech.X_from_kesten_f(0.0) == pytest.approx(0.25)
    assert mech.X_from_kesten_f(1.0) == pytest.approx(1.0)


def test_kesten_f_conversion_roundtrips(mech):
    f = np.array([0.0, 0.3, 0.645, 1.0])
    assert np.allclose(mech.kesten_f_from_X(mech.X_from_kesten_f(f)), f)


#: Two stations from Kesten's own program output
#: (docs/kesten_claude/vapor_reference.csv), feed at TF = 530 degR.
KESTEN_STATIONS = [
    # (z_ft,      T_degR,      f_kesten)
    (0.17155202, 1950.7243, 0.59806513),
    (0.24999999, 1905.3842, 0.64541664),
]


@pytest.mark.parametrize("z_ft,T_degR,f", KESTEN_STATIONS)
def test_adiabatic_temperature_matches_kesten_reactor_output(mech, z_ft, T_degR, f):
    """Agree with Kesten's 1967 code to a few percent at his own conditions.

    Independent thermodynamics (his tabulated cp, our NASA-7 polynomials) and
    an independent implementation, so a few tens of K is a good result. The
    conversion from his dissociation fraction is essential -- without it the
    error is around 90 K rather than 20.
    """
    T_ref = T_degR * RANKINE_TO_KELVIN
    T_feed = 530.0 * RANKINE_TO_KELVIN
    T_mine = mech.adiabatic_temperature(mech.X_from_kesten_f(f), T_feed=T_feed)
    assert T_mine == pytest.approx(T_ref, rel=0.03)


def test_skipping_the_f_conversion_is_visibly_wrong(mech):
    """Guard the trap itself: using f as X must disagree badly."""
    T_feed = 530.0 * RANKINE_TO_KELVIN
    T_ref = 1905.3842 * RANKINE_TO_KELVIN
    T_naive = mech.adiabatic_temperature(0.64541664, T_feed=T_feed)
    assert abs(T_naive - T_ref) > 50.0


# -- rate laws --------------------------------------------------------------
def test_catalytic_rate_increases_with_solid_temperature(mech):
    rate = CatalyticRate(1e10, 100e6, reactant="NH3")
    Y = mech.composition_at_dissociation(0.3)
    assert rate(_state(mech, Y, T_solid=1200.0)) > rate(_state(mech, Y, T_solid=800.0))


def test_homogeneous_rate_follows_gas_temperature(mech):
    """Unlike the catalytic paths, which follow the solid."""
    rate = HomogeneousRate(1e10, 100e6, reactant="N2H4")
    Y = mech.inlet_composition()
    hot_gas = rate(_state(mech, Y, T_gas=1200.0, T_solid=400.0))
    cold_gas = rate(_state(mech, Y, T_gas=400.0, T_solid=1200.0))
    assert hot_gas > cold_gas


def test_hydrogen_inhibits_ammonia_dissociation(mech):
    """The self-limiting feedback: more H2 must slow the reaction."""
    rate = mech.reactions[1].rate
    base = mech.composition_at_dissociation(0.3)
    lo = rate(_state(mech, base))
    more_h2 = base.copy()
    i_h2, i_n2 = mech.mixture.index("H2"), mech.mixture.index("N2")
    more_h2[i_h2] *= 3.0
    more_h2[i_n2] -= base[i_h2] * 2.0
    assert rate(_state(mech, more_h2 / more_h2.sum())) < lo


def test_hydrogen_floor_keeps_the_rate_finite(mech):
    """1/C_H2^1.6 diverges at zero hydrogen; the floor must contain it."""
    Y = np.zeros(mech.mixture.n)
    Y[mech.mixture.index("NH3")] = 0.5
    Y[mech.mixture.index("N2")] = 0.5          # no H2 at all
    r = mech.rates(T_gas=900.0, T_solid=1100.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    assert np.all(np.isfinite(r))


def test_series_resistance_is_bounded_by_both_branches():
    a = np.array([1.0, 100.0, 1e6, 0.0])
    b = np.array([100.0, 1.0, 1.0, 5.0])
    s = series(a, b)
    assert np.all(s <= a + 1e-12) and np.all(s <= b + 1e-12)


def test_kesten_mass_transfer_rises_with_mass_flux(mech):
    mt = KestenMassTransfer(0.17e-3 * 0.09290304, "NH3")
    Y = mech.composition_at_dissociation(0.3)
    assert mt(_state(mech, Y, G=50.0)) > mt(_state(mech, Y, G=5.0))


def test_diffusivity_scales_with_temperature_and_pressure():
    mt = KestenMassTransfer(1e-5, "NH3")
    assert mt.diffusivity(mt.T_REF, mt.P_REF) == pytest.approx(1e-5, rel=1e-12)
    assert mt.diffusivity(2 * mt.T_REF, mt.P_REF) == pytest.approx(1e-5 * 2 ** 1.823, rel=1e-9)
    assert mt.diffusivity(mt.T_REF, 2 * mt.P_REF) == pytest.approx(0.5e-5, rel=1e-12)


def test_rates_are_nonnegative_and_finite(mech):
    r = mech.rates(T_gas=600.0, T_solid=900.0, p=8e5,
                   Y=mech.inlet_composition(), G=20.0, bed=FakeBed())
    assert np.all(np.isfinite(r)) and np.all(r >= 0.0)


def test_no_reaction_without_reactant(mech):
    Y = mech.composition_at_dissociation(1.0)
    r = mech.rates(T_gas=1000.0, T_solid=1000.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    assert np.allclose(r, 0.0, atol=1e-12)


def test_hydrazine_step_is_diffusion_limited(mech):
    """Kesten's claim, and the reason his arbitrary Ea1 does not matter.

    If raising the pre-exponential moves the rate, kinetics is controlling and
    the least trustworthy number in the model is setting the answer.
    """
    Y = mech.inlet_composition()
    kw = dict(T_gas=800.0, T_solid=1100.0, p=8e5, Y=Y, G=20.0, bed=FakeBed())
    slow = HydrazineShell405(A_N2H4_cat=1.0e10).rates(**kw)[0]
    fast = HydrazineShell405(A_N2H4_cat=1.0e13).rates(**kw)[0]
    assert fast == pytest.approx(slow, rel=0.02)


def test_homogeneous_path_can_be_disabled(mech):
    assert len(HydrazineShell405(include_homogeneous=False).reactions) == 2
    assert len(mech.reactions) == 3


def test_inhibition_order_is_selectable(mech):
    """Melton reports 1.0, Logan and Kemball 1.6, and Kesten's own Fortran
    disagrees with itself. Both must be reachable."""
    for n in (1.0, 1.6):
        m = HydrazineShell405(nH2=n)
        assert m.reactions[1].rate.inhibitor_order == pytest.approx(n)


# -- registry ---------------------------------------------------------------
def test_registry_returns_a_mechanism():
    assert isinstance(get_mechanism("hydrazine"), HydrazineShell405)


def test_registry_rejects_unknown_names():
    with pytest.raises(ValueError, match="unknown mechanism"):
        get_mechanism("unobtainium")
