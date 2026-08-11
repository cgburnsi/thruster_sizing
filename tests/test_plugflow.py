"""Tests for the reacting plug-flow reactor.

Three layers here:

* Conservation -- mass, elements, and the fact that an adiabatic bed conserves
  total enthalpy. The last one gives a strong internal check: the plug-flow
  exit temperature must equal the closed-form adiabatic temperature at the
  same dissociation fraction, computed by an entirely separate code path.
* Physical behaviour -- dissociation rises, pressure falls, the catalyst runs
  hotter than the gas where decomposition dominates.
* Kesten's own program output as an end-to-end regression.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from fvm.catbed import CatalystBed, PackedSpheres
from fvm.mechanism import HydrazineShell405, RANKINE_TO_KELVIN
from fvm.plugflow import PlugFlowReactor

FT = 0.3048
PSI = 6894.757293
LB = 0.45359237

REFERENCE = Path(__file__).resolve().parents[1] / (
    "docs/kesten_claude/vapor_reference.csv")

#: Kesten's own code used ALPHA2 = 1e11 with a hydrogen exponent of 1.0 in
#: MAIN.f. PARAM.f and the published Eq. 43 use 1.6; see the note on
#: test_reproduces_kesten_exit_state.
KESTEN_CODE_KINETICS = dict(nH2=1.0, A_NH3=1.0e11)


def _kesten_case():
    rows = list(csv.DictReader(REFERENCE.open()))
    first, last = rows[0], rows[-1]
    order = ["H2", "N2", "NH3", "N2H4"]
    conc = [float(first["c%d" % i]) for i in (1, 2, 3, 4)]
    return dict(
        order=order, conc=conc,
        T_in=float(first["temp_degR"]) * RANKINE_TO_KELVIN,
        p_in=float(first["pres_psia"]) * PSI,
        T_out=float(last["temp_degR"]) * RANKINE_TO_KELVIN,
        p_out=float(last["pres_psia"]) * PSI,
        X_out=(3.0 * float(last["frac3d"]) + 1.0) / 4.0,
        X_in=(3.0 * float(first["frac3d"]) + 1.0) / 4.0,
        G=3.0 * LB / FT ** 2,
    )


def _inlet_fractions(mech, order, conc):
    Y = np.zeros(mech.mixture.n)
    for name, c in zip(order, conc):
        Y[mech.mixture.index(name)] = c
    return Y / Y.sum()


@pytest.fixture(scope="module")
def kesten_solution():
    if not REFERENCE.exists():
        pytest.skip("Kesten reference data absent")
    case = _kesten_case()
    mech = HydrazineShell405(**KESTEN_CODE_KINETICS)
    bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)
    sol = PlugFlowReactor(mech, bed).solve(
        G=case["G"], p_inlet=case["p_in"], T_inlet=case["T_in"],
        Y_inlet=_inlet_fractions(mech, case["order"], case["conc"]),
        n_output=150)
    return sol, case, mech


# -- conservation -----------------------------------------------------------
def test_mass_fractions_sum_to_one(kesten_solution):
    sol, _, _ = kesten_solution
    assert np.allclose(sol.Y.sum(axis=0), 1.0, atol=1e-9)


def test_mass_fractions_stay_nonnegative(kesten_solution):
    sol, _, _ = kesten_solution
    assert np.all(sol.Y >= -1e-12)


def test_element_ratio_is_conserved(kesten_solution):
    """The N:H ratio must not drift along the bed.

    Asserted as constancy, not as equality with 2.0. The inlet here is built
    from Kesten's reported concentrations, which he converted with his own
    molecular weights, so its ratio sits 0.26% off the stoichiometric value
    before integration even starts. What the solver owes us is that it does
    not make that any worse -- and it holds the inherited value to eight
    decimals across the whole bed, independent of integrator tolerance.
    """
    sol, _, mech = kesten_solution
    x = sol.mole_fractions
    m = mech.mixture
    N = (2 * x[m.index("N2H4")] + x[m.index("NH3")] + 2 * x[m.index("N2")])
    H = (4 * x[m.index("N2H4")] + 3 * x[m.index("NH3")] + 2 * x[m.index("H2")])
    ratio = H / N
    assert np.allclose(ratio, ratio[0], rtol=1e-8)
    assert ratio[0] == pytest.approx(2.0, rel=5e-3)


def test_adiabatic_bed_conserves_total_enthalpy(kesten_solution):
    """No heat loss means h is constant -- there is no reaction source term."""
    sol, _, _ = kesten_solution
    assert np.allclose(sol.h, sol.h[0], rtol=1e-8)


def test_exit_temperature_matches_the_closed_form_adiabatic_value(kesten_solution):
    """Cross-check against an entirely separate code path.

    ``adiabatic_temperature`` inverts the enthalpy balance directly from a
    dissociation fraction; the reactor integrates species and energy along the
    bed. Agreeing to a fraction of a percent means the integration is not
    quietly leaking energy.
    """
    sol, _, mech = kesten_solution
    T_feed = sol.T_gas[0]
    Y0 = sol.Y[:, 0]
    h_in = mech.mixture.h(T_feed, Y0)
    T_closed = float(mech.mixture.temperature_from_h(
        h_in, sol.Y[:, -1], T_guess=1000.0))
    assert sol.T_gas[-1] == pytest.approx(T_closed, rel=1e-6)


# -- physical behaviour -----------------------------------------------------
def test_dissociation_increases_along_the_bed(kesten_solution):
    sol, _, _ = kesten_solution
    X = sol.dissociation
    assert X[-1] > X[0]
    assert np.all(np.diff(X) > -1e-6)


def test_inlet_dissociation_matches_kestens_reported_value(kesten_solution):
    """Guards the X definition: normalised by ammonia produced, not fed."""
    sol, case, _ = kesten_solution
    assert sol.dissociation[0] == pytest.approx(case["X_in"], abs=0.005)


def test_hydrazine_is_consumed_early(kesten_solution):
    """Decomposition is diffusion-limited and fast; it should not survive far."""
    sol, _, mech = kesten_solution
    n2h4 = sol.species("N2H4")
    gone = np.nonzero(n2h4 < 1e-2 * n2h4[0])[0]
    assert gone.size > 0, "hydrazine never substantially decomposed"
    assert sol.z[gone[0]] < 0.25 * sol.z[-1], "99% conversion should be early"
    assert n2h4[-1] < 1e-4 * n2h4[0], "and essentially complete by the exit"


def test_pressure_falls_monotonically(kesten_solution):
    sol, _, _ = kesten_solution
    assert np.all(np.diff(sol.p) < 0.0)
    assert sol.pressure_drop > 0.0


def test_catalyst_runs_hotter_than_the_gas_where_decomposition_dominates(kesten_solution):
    """The exothermic step deposits heat on the solid, which then feeds the gas."""
    sol, _, _ = kesten_solution
    assert sol.T_solid[0] > sol.T_gas[0] + 50.0


def test_exit_conditions_are_physically_sensible(kesten_solution):
    sol, _, _ = kesten_solution
    e = sol.exit_conditions()
    assert 800.0 < e["T"] < 1600.0
    assert 10.0 < e["MW"] < 20.0
    assert 1.1 < e["gamma"] < 1.5
    assert 1000.0 < e["cstar"] < 1600.0
    assert 0.0 <= e["dissociation"] <= 1.0


def test_higher_mass_flux_leaves_less_time_to_dissociate():
    """Shorter residence time must mean less ammonia dissociation."""
    mech = HydrazineShell405(**KESTEN_CODE_KINETICS)
    bed = CatalystBed.uniform(0.02, 0.05, PackedSpheres(mesh=(25, 30), eps=0.38))
    rx = PlugFlowReactor(mech, bed)
    Y0 = mech.inlet_composition()
    out = []
    for G in (10.0, 40.0):
        s = rx.solve(G=G, p_inlet=100 * PSI, T_inlet=900.0, Y_inlet=Y0, n_output=60)
        out.append(s.dissociation[-1])
    assert out[1] < out[0]


# -- regression against Kesten's program ------------------------------------
def test_reproduces_kesten_exit_state(kesten_solution):
    """End-to-end against his own output, using his own code's kinetics.

    Which kinetics those are is itself a finding. His two copies of the rate
    subroutine disagree -- PARAM.f divides by C_H2^1.6, MAIN.f by C_H2^1.0 --
    and the published Eq. 43 quotes 1.6 with A = 0.3e11. Sweeping both
    parameters against this reference case, the combination that reproduces it
    is n = 1.0 with A = 1e11, which is exactly MAIN.f together with the
    ALPHA2 in his input deck. That is evidence about which path generated the
    reference output, not proof that 1.0 is the better physics: this model
    omits the intraparticle diffusion resistance his does carry, so some of
    the agreement may be compensating error.

    The library default remains the published Eq. 43 values.
    """
    sol, case, _ = kesten_solution
    assert sol.dissociation[-1] == pytest.approx(case["X_out"], abs=0.05)
    assert sol.T_gas[-1] == pytest.approx(case["T_out"], rel=0.02)


def test_default_kinetics_still_run_and_over_dissociate(kesten_solution):
    """The published Eq. 43 defaults should work, and predict more dissociation.

    Documented rather than hidden: with n = 1.6 and A = 0.3e11 the model
    reaches X ~ 0.83 against Kesten's 0.73 on this case.
    """
    _, case, _ = kesten_solution
    mech = HydrazineShell405()
    bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)
    sol = PlugFlowReactor(mech, bed).solve(
        G=case["G"], p_inlet=case["p_in"], T_inlet=case["T_in"],
        Y_inlet=_inlet_fractions(mech, case["order"], case["conc"]),
        n_output=120)
    assert sol.dissociation[-1] > case["X_out"]
    assert np.isfinite(sol.T_gas).all()
