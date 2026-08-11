"""Tests for the system-level bed/nozzle coupling."""
import numpy as np
import pytest

from fvm.catbed import CatalystBed, PackedSpheres
from fvm.geometry import ConicalNozzle
from fvm.mechanism import HydrazineShell405
from fvm.thruster import ThrusterSystem, perfect_gas_from_bed, vapor_region_inlet

MDOT = 0.0443e-3          # kg/s, roughly 100 mN of hydrazine

#: These are system-level integration tests -- each solve runs the bed
#: integration to convergence -- so the whole module is slow by construction.
pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def mech():
    return HydrazineShell405(nH2=1.0, A_NH3=1.0e11)


@pytest.fixture(scope="module")
def system(mech):
    bed = CatalystBed.uniform(8.3e-3, 20e-3,
                              PackedSpheres(mesh=(25, 30), eps=0.38))
    nozzle = ConicalNozzle(4.15e-3, 0.145e-3, 100.0, 4e-3)
    return ThrusterSystem(mech, bed, nozzle, p_ambient=666.6, T_vapor=455.6)


@pytest.fixture(scope="module")
def solution(system):
    return system.solve_for_mdot(MDOT)


# -- the vapour-region inlet ------------------------------------------------
def test_vapor_inlet_conserves_enthalpy_from_the_liquid_feed(mech):
    """The defining property: h(T_vapor, Y) must equal the liquid feed's h."""
    Y, f = vapor_region_inlet(mech, 455.6, 298.15)
    assert mech.mixture.h(455.6, Y) == pytest.approx(
        mech.inlet_enthalpy(298.15), rel=1e-9)
    assert 0.0 < f < 1.0


def test_vapor_inlet_extent_is_close_to_kestens(mech):
    """Kesten's own vapour region begins ~38% decomposed at 455.6 K.

    Nothing here is fitted to that -- the extent falls out of requiring energy
    conservation through vaporisation -- so landing within a few points of his
    is an independent check on the thermodynamics.
    """
    _, f = vapor_region_inlet(mech, 455.6, 298.15)
    assert f == pytest.approx(0.38, abs=0.06)


def test_hotter_vapor_inlet_needs_more_predecomposition(mech):
    _, cold = vapor_region_inlet(mech, 420.0, 298.15)
    _, hot = vapor_region_inlet(mech, 500.0, 298.15)
    assert hot > cold


def test_vapor_inlet_rejects_unreachable_temperatures(mech):
    """Only the upper bound is reachable for hydrazine.

    The lower guard would need undecomposed vapour to carry *less* enthalpy
    than the liquid, which for hydrazine cannot happen at any physical
    temperature: it is an endothermic compound, so vaporisation and formation
    both push the vapour's enthalpy above the liquid's. The guard stays for
    propellants where that is not true.

    The upper bound is the adiabatic temperature at zero dissociation, 1646 K
    -- decomposing every last molecule cannot get the vapour hotter than that.
    """
    with pytest.raises(ValueError, match="too high"):
        vapor_region_inlet(mech, 1800.0, 298.15)


def test_pure_vapour_inlet_would_overstate_the_flame_temperature(mech):
    """Why the energy-consistent inlet matters, pinned as a test.

    Feeding pure hydrazine vapour at T_vapor hands the bed its vaporisation
    enthalpy for free. At complete dissociation that puts the adiabatic
    temperature hundreds of kelvin above the correct value.
    """
    Y_pure = mech.inlet_composition()
    Y_real, _ = vapor_region_inlet(mech, 455.6, 298.15)
    products = mech.composition_at_dissociation(1.0)
    T_pure = float(mech.mixture.temperature_from_h(
        mech.mixture.h(455.6, Y_pure), products, T_guess=1200.0))
    T_real = float(mech.mixture.temperature_from_h(
        mech.mixture.h(455.6, Y_real), products, T_guess=900.0))
    assert T_pure > T_real + 300.0
    assert T_real == pytest.approx(mech.adiabatic_temperature(1.0), rel=1e-6)


# -- gas handoff ------------------------------------------------------------
def test_perfect_gas_matches_the_bed_exit(solution, mech):
    e = solution.chamber
    gas = perfect_gas_from_bed(e, mech)
    assert gas.gamma == pytest.approx(e["gamma"])
    assert gas.MW == pytest.approx(e["MW"])
    assert gas.R == pytest.approx(e["R"], rel=1e-9)


# -- the coupled solution ---------------------------------------------------
def test_pressures_are_ordered(solution):
    assert solution.p_feed > solution.p_chamber > solution.p_ambient
    assert solution.bed_pressure_drop > 0.0


def test_chamber_temperature_matches_the_closed_form(solution, mech):
    """The reactor and the closed-form enthalpy balance must agree exactly."""
    T_closed = mech.adiabatic_temperature(
        solution.chamber["dissociation"], T_feed=298.15)
    assert solution.chamber["T"] == pytest.approx(T_closed, rel=1e-4)


def test_throat_passes_the_requested_mass_flow(solution):
    """The whole point of the iteration: the choked throat must be consistent."""
    e = solution.chamber
    mdot_throat = solution.Cd * e["p"] * np.pi * \
        solution.contour.r_throat ** 2 / e["cstar"]
    assert mdot_throat == pytest.approx(solution.mdot, rel=1e-4)


def test_performance_is_in_the_hydrazine_band(solution):
    """Monopropellant hydrazine runs about 220-235 s."""
    assert 200.0 < solution.Isp < 245.0
    assert 1150.0 < solution.chamber["cstar"] < 1400.0
    assert 0.0 < solution.chamber["dissociation"] <= 1.0


def test_solution_round_trips(system, solution):
    """Solve for feed pressure, then invert it and recover the mass flow."""
    back = system.solve_for_feed_pressure(solution.p_feed)
    assert back.mdot == pytest.approx(solution.mdot, rel=1e-4)


def test_more_flow_needs_more_feed_pressure(system, solution):
    higher = system.solve_for_mdot(MDOT * 1.5)
    assert higher.p_feed > solution.p_feed
    assert higher.thrust > solution.thrust


def test_discharge_coefficient_raises_the_required_chamber_pressure(mech, solution):
    """A real throat passes less, so it needs more pressure for the same flow.

    This is the feedback path from the CFD solver: run it, get Cd, re-converge.
    """
    bed = CatalystBed.uniform(8.3e-3, 20e-3,
                              PackedSpheres(mesh=(25, 30), eps=0.38))
    nozzle = ConicalNozzle(4.15e-3, 0.145e-3, 100.0, 4e-3)
    viscous = ThrusterSystem(mech, bed, nozzle, p_ambient=666.6,
                             T_vapor=455.6, Cd=0.93).solve_for_mdot(MDOT)
    assert viscous.p_chamber > solution.p_chamber
    assert viscous.p_chamber / solution.p_chamber == pytest.approx(1 / 0.93, rel=0.02)


def test_report_renders(solution):
    text = solution.report(quiet=True)
    assert "Chamber pressure" in text and "Specific impulse" in text
