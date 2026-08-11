"""Tests for multi-species thermodynamics.

The point of this file is to catch transcription errors in the NASA-7
coefficient table. A wrong digit shifts cp by a few percent and nothing
crashes, so every built-in species is checked against independently known
property values rather than against itself.

Reference values are standard textbook/JANAF numbers, quoted to the precision
they are usually tabulated at, with tolerances set accordingly.
"""
import numpy as np
import pytest

from fvm import chem


# -- cp against known values [J/(kg*K)] -------------------------------------
CP_REFERENCE = {
    #  species: [(T, cp [J/(kg*K)], rel_tol), ...]
    "N2":  [(300.0, 1040.0, 0.01), (1000.0, 1167.0, 0.01), (2000.0, 1284.0, 0.015)],
    "H2":  [(300.0, 14310.0, 0.01), (1000.0, 14980.0, 0.01)],
    # JANAF molar values converted at MW 17.0306: Cp = 35.65 J/(mol*K) at
    # 298 K and 56.3 at 1000 K. An earlier revision of this test quoted 2190
    # and 3000 from memory and failed -- the table was right and the test was
    # wrong, which is worth remembering before "fixing" a coefficient.
    "NH3": [(298.15, 2093.0, 0.01), (1000.0, 3305.0, 0.01)],
    "H2O": [(300.0, 1864.0, 0.01), (1000.0, 2290.0, 0.02)],
    "O2":  [(300.0, 918.0, 0.01), (1000.0, 1090.0, 0.02)],
}


@pytest.mark.parametrize("name", sorted(CP_REFERENCE))
def test_cp_matches_reference_values(name):
    s = chem.species(name)
    for T, cp_ref, rtol in CP_REFERENCE[name]:
        assert s.cp(T) == pytest.approx(cp_ref, rel=rtol), \
            f"{name} cp at {T} K is {s.cp(T):.1f}, expected ~{cp_ref}"


# -- enthalpy of formation at 298.15 K [J/kmol] -----------------------------
HF_REFERENCE = {
    "N2": (0.0, 2.0e6),            # elements in standard state are exactly zero
    "H2": (0.0, 2.0e6),
    "O2": (0.0, 2.0e6),
    "NH3": (-45.94e6, 1.0e6),      # -45.94 kJ/mol
    "H2O": (-241.83e6, 1.0e6),     # -241.83 kJ/mol (gas)
    "N2H4": (95.35e6, 3.0e6),      # +95.35 kJ/mol (gas) -- endothermic compound
}


@pytest.mark.parametrize("name", sorted(HF_REFERENCE))
def test_enthalpy_of_formation(name):
    ref, atol = HF_REFERENCE[name]
    hf = chem.species(name).h_formation()
    assert hf == pytest.approx(ref, abs=atol), \
        f"{name} h_f = {hf / 1e6:.2f} MJ/kmol, expected ~{ref / 1e6:.2f}"


def test_hydrazine_is_endothermic():
    """The whole monopropellant works because N2H4 stores energy.

    If this flips sign, the bed model will predict cooling instead of a
    temperature rise, so it is worth asserting on its own.
    """
    assert chem.species("N2H4").h_formation() > 0.0


def test_molecular_weights():
    for name, mw in (("N2", 28.0134), ("H2", 2.01594), ("NH3", 17.0306),
                     ("H2O", 18.0153), ("O2", 31.9988), ("N2H4", 32.0452)):
        assert chem.species(name).MW == pytest.approx(mw, rel=1e-4)


def test_gamma_of_diatomics_near_room_temperature():
    """N2 and H2 should both sit close to the diatomic value of 7/5."""
    for name in ("N2", "H2"):
        m = chem.mixture([name], [1.0])
        assert m.gamma(300.0) == pytest.approx(1.40, abs=0.01)


#: Published NASA-7 fits are normally continuous at T_mid to round-off; the
#: NH3 entry carries a genuine 0.39% cp mismatch that comes with the data. The
#: 1% gate tolerates that while still catching a mistyped exponent, which
#: shifts cp by several percent or more.
_CONTINUITY_TOL = 0.01


def test_cp_is_continuous_across_the_polynomial_switch():
    """A mismatched low/high coefficient set shows up as a jump at T_mid."""
    for name in sorted(chem._DATA):
        s = chem.species(name)
        Tm = s.T_mid
        lo, hi = s.cp(Tm - 1e-3), s.cp(Tm + 1e-3)
        assert abs(hi - lo) / lo < _CONTINUITY_TOL, \
            f"{name} cp jumps at T_mid: {lo:.2f} -> {hi:.2f}"


def test_enthalpy_is_continuous_across_the_polynomial_switch():
    """Compared against cp*T, which is the scale an enthalpy error matters on.

    A relative test on h itself is meaningless: h carries the formation
    enthalpy, so it passes through zero for the elements.
    """
    for name in sorted(chem._DATA):
        s = chem.species(name)
        Tm = s.T_mid
        jump = abs(s.h(Tm + 1e-3) - s.h(Tm - 1e-3))
        assert jump < _CONTINUITY_TOL * s.cp(Tm) * Tm, \
            f"{name} h jumps at T_mid by {jump:.1f} J/kg"


# -- mixture behaviour ------------------------------------------------------
def test_mixture_molecular_weight():
    m = chem.mixture(["N2", "H2"])
    Y = m.from_moles({"N2": 1.0, "H2": 3.0})
    expected = (28.0134 + 3 * 2.01594) / 4.0
    assert m.MW(Y) == pytest.approx(expected, rel=1e-10)


def test_mole_and_mass_fraction_roundtrip():
    m = chem.mixture(["N2", "H2", "NH3"])
    moles = {"N2": 0.5, "H2": 1.5, "NH3": 2.0}
    Y = m.from_moles(moles)
    x = m.mole_fractions(Y)
    tot = sum(moles.values())
    for name, n in moles.items():
        assert x[m.index(name)] == pytest.approx(n / tot, rel=1e-10)


def test_mixture_cp_is_mass_weighted():
    m = chem.mixture(["N2", "H2"])
    Y = np.array([0.75, 0.25])
    T = 800.0
    expected = 0.75 * chem.species("N2").cp(T) + 0.25 * chem.species("H2").cp(T)
    assert m.cp(T, Y) == pytest.approx(expected, rel=1e-12)


def test_temperature_from_enthalpy_inverts_h():
    m = chem.mixture(["N2", "H2", "NH3"])
    Y = m.from_moles({"N2": 1.0, "H2": 2.0, "NH3": 1.0})
    for T_true in (400.0, 900.0, 1500.0, 2500.0):
        h = m.h(T_true, Y)
        assert m.temperature_from_h(h, Y, T_guess=300.0) == pytest.approx(T_true, rel=1e-6)


def test_temperature_from_enthalpy_is_vectorised():
    m = chem.mixture(["N2", "H2"])
    Y = m.from_moles({"N2": 1.0, "H2": 1.0})
    T_true = np.array([500.0, 1200.0, 2200.0])
    T = m.temperature_from_h(m.h(T_true, Y), Y, T_guess=1000.0)
    assert np.allclose(T, T_true, rtol=1e-6)


def test_gamma_falls_as_temperature_rises():
    """Vibrational modes activate, so gamma must decrease monotonically."""
    m = chem.mixture(["N2"], [1.0])
    g = np.array([m.gamma(T) for T in (300.0, 800.0, 1500.0, 3000.0)])
    assert np.all(np.diff(g) < 0.0)


def test_unknown_species_raises_with_a_useful_message():
    with pytest.raises(KeyError, match="no built-in data"):
        chem.species("Unobtainium")
