"""Tests for catalyst bed geometry and pressure drop.

The external check is Kesten's own pressure profile: his program reports the
bed dropping 100 to 71.68 psia, and his reference output carries enough
composition data to reconstruct the density along the bed. Ergun should
reproduce that to within its documented accuracy, which is not tight -- the
correlation is generally quoted at +/-25% and worse for non-spherical packing.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from fvm.catbed import (CatalystBed, PackedCylinders, PackedSpheres,
                        ReticulatedFoam, SIEVE_OPENING, mesh_to_diameter)
from fvm.mechanism import (HydrazineShell405, LB_FT3_TO_KG_M3,
                           RANKINE_TO_KELVIN)

FT = 0.3048
PSI = 6894.757293
LB = 0.45359237


# -- sieve sizes ------------------------------------------------------------
def test_single_mesh_returns_its_opening():
    assert mesh_to_diameter(30) == pytest.approx(600e-6)


def test_mesh_cut_is_the_geometric_mean():
    """25-30 mesh is material between 0.600 and 0.710 mm."""
    d = mesh_to_diameter(25, 30)
    assert d == pytest.approx(np.sqrt(710e-6 * 600e-6))
    assert 600e-6 < d < 710e-6


def test_sieve_openings_decrease_with_mesh_number():
    n = sorted(SIEVE_OPENING)
    assert np.all(np.diff([SIEVE_OPENING[k] for k in n]) < 0)


def test_reversed_mesh_cut_is_rejected():
    with pytest.raises(ValueError, match="finer mesh"):
        mesh_to_diameter(30, 25)


def test_unknown_mesh_gives_a_useful_error():
    with pytest.raises(KeyError, match="no sieve opening"):
        mesh_to_diameter(999)


# -- packing geometry -------------------------------------------------------
def test_sphere_sphericity_is_unity_and_area_matches_the_identity():
    m = PackedSpheres(diameter=1e-3, eps=0.4)
    assert m.sphericity == pytest.approx(1.0)
    assert m.d_p == pytest.approx(1e-3)
    assert m.a_v == pytest.approx(6.0 * 0.6 / 1e-3)


def test_equal_cylinder_sphericity_is_the_textbook_value():
    """L = D cylinders have sphericity 0.874, computed here, not tabulated."""
    m = PackedCylinders(0.125 * 0.0254)
    assert m.sphericity == pytest.approx(0.874, abs=0.003)


def test_cylinder_volume_equivalent_diameter():
    d = 2e-3
    m = PackedCylinders(d, d)
    volume = np.pi * d ** 3 / 4.0
    assert m.d_v == pytest.approx((6 * volume / np.pi) ** (1 / 3))


def test_finer_packing_gives_more_surface_area():
    coarse = PackedSpheres(diameter=2e-3, eps=0.4)
    fine = PackedSpheres(diameter=0.5e-3, eps=0.4)
    assert fine.a_v == pytest.approx(4 * coarse.a_v)


def test_foam_surface_area_matches_the_thruster_sizing_correlation():
    """a_v = 1300 (PPI/20)^1.3, the expression already in thruster_sizing.py."""
    f = ReticulatedFoam(ppi=60)
    assert f.a_v == pytest.approx(1300.0 * 3.0 ** 1.3, rel=1e-9)
    assert f.a_v == pytest.approx(5422.0, rel=0.01)


def test_foam_is_mostly_void():
    assert ReticulatedFoam(ppi=60).eps > 0.8


# -- Ergun ------------------------------------------------------------------
def test_pressure_gradient_is_positive_and_grows_with_flux():
    m = PackedSpheres(diameter=1e-3, eps=0.4)
    lo = m.pressure_gradient(G=1.0, rho=1.0, mu=3e-5)
    hi = m.pressure_gradient(G=10.0, rho=1.0, mu=3e-5)
    assert 0 < lo < hi


def test_viscous_term_dominates_at_low_reynolds_number():
    """Creeping flow must be linear in G."""
    m = PackedSpheres(diameter=1e-4, eps=0.4)
    kw = dict(rho=1.0, mu=1e-3)
    assert m.reynolds(G=1e-3, mu=1e-3) < 1.0
    r = m.pressure_gradient(G=2e-3, **kw) / m.pressure_gradient(G=1e-3, **kw)
    assert r == pytest.approx(2.0, rel=0.02)


def test_inertial_term_dominates_at_high_reynolds_number():
    """Fully inertial flow must go as G^2."""
    m = PackedSpheres(diameter=3e-3, eps=0.4)
    kw = dict(rho=1.0, mu=1e-5)
    assert m.reynolds(G=100.0, mu=1e-5) > 1e4
    r = m.pressure_gradient(G=200.0, **kw) / m.pressure_gradient(G=100.0, **kw)
    assert r == pytest.approx(4.0, rel=0.02)


def test_pressure_drop_rises_steeply_as_voidage_falls():
    """Ergun goes as (1-eps)/eps^3, so the bed tightens fast."""
    kw = dict(G=10.0, rho=1.0, mu=3e-5)
    loose = PackedSpheres(diameter=1e-3, eps=0.45).pressure_gradient(**kw)
    tight = PackedSpheres(diameter=1e-3, eps=0.35).pressure_gradient(**kw)
    assert tight > 2.0 * loose


# -- bed assembly -----------------------------------------------------------
def test_bed_geometry():
    bed = CatalystBed.uniform(0.02, 0.05, PackedSpheres(diameter=1e-3))
    assert bed.length == pytest.approx(0.05)
    assert bed.area == pytest.approx(np.pi * 0.01 ** 2)
    assert bed.volume == pytest.approx(bed.area * 0.05)


def test_sections_are_selected_by_position():
    a = PackedSpheres(diameter=1e-3, name="fine")
    b = PackedCylinders(3e-3, name="coarse")
    bed = CatalystBed(0.02, [(0.01, a), (0.04, b)])
    assert bed.medium_at(0.0).name == "fine"
    assert bed.medium_at(0.005).name == "fine"
    assert bed.medium_at(0.02).name == "coarse"
    assert bed.medium_at(bed.length).name == "coarse"


def test_empty_or_negative_sections_are_rejected():
    with pytest.raises(ValueError):
        CatalystBed(0.02, [])
    with pytest.raises(ValueError):
        CatalystBed(0.02, [(-0.01, PackedSpheres(diameter=1e-3))])


def test_kesten_standard_bed_matches_his_description():
    """25-30 mesh for the first 0.2 in, 1/8 x 1/8 in pellets after."""
    bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)
    assert bed.length == pytest.approx(0.25 * FT)
    assert bed.sections[0][0] == pytest.approx(0.2 * 0.0254)
    assert bed.medium_at(0.001).d_p == pytest.approx(mesh_to_diameter(25, 30) * 1.0)
    assert bed.medium_at(0.05).d_p > bed.medium_at(0.001).d_p


def test_kesten_bed_rejects_a_length_shorter_than_its_entry_zone():
    with pytest.raises(ValueError, match="shorter than"):
        CatalystBed.kesten_standard(diameter=0.02, length=0.1 * 0.0254)


def test_mass_flux_and_loading_units():
    bed = CatalystBed.uniform(0.02, 0.05, PackedSpheres(diameter=1e-3))
    mdot = 0.01
    assert bed.mass_flux(mdot) == pytest.approx(mdot / bed.area)
    # 1 kg/m^2-s is 0.00142233 lb/in^2-s
    assert bed.loading_lb_in2_s(mdot) == pytest.approx(
        bed.mass_flux(mdot) * 0.00142233)


def test_pressure_drop_accepts_callable_properties():
    """Density varies by orders of magnitude in a reacting bed."""
    bed = CatalystBed.uniform(0.02, 0.05, PackedSpheres(diameter=1e-3))
    const = bed.pressure_drop(G=10.0, rho=1.0, mu=3e-5)
    varying = bed.pressure_drop(G=10.0, rho=lambda z: 1.0, mu=lambda z: 3e-5)
    assert const == pytest.approx(varying, rel=1e-6)


# -- validation against Kesten's pressure profile ---------------------------
REFERENCE = Path(__file__).resolve().parents[1] / (
    "docs/kesten_claude/vapor_reference.csv")


def _kesten_profile():
    rows = list(csv.DictReader(REFERENCE.open()))
    d = np.array([[float(r["z_ft"]), float(r["temp_degR"]), float(r["pres_psia"])]
                  + [float(r["c%d" % k]) for k in (1, 2, 3, 4)] for r in rows])
    d = d[np.argsort(d[:, 0])]
    return (d[:, 0] * FT, d[:, 1] * RANKINE_TO_KELVIN, d[:, 2] * PSI,
            d[:, 3:].sum(axis=1) * LB_FT3_TO_KG_M3)


@pytest.mark.skipif(not REFERENCE.exists(), reason="Kesten reference data absent")
def test_ergun_reproduces_kesten_pressure_drop_on_the_resolved_segment():
    """Check the one segment where his output actually resolves the bed.

    His six stations put four inside the first 0.4 mm of a 76 mm bed, so the
    density profile over most of the length is unresolved and interpolating it
    is guesswork. The final 24 mm is the exception: real length, and density
    varying by only 9%. Ergun should land within its stated accuracy there.

    It comes out around 0.8 of measured. The shortfall is consistent across
    every segment despite a 6x density range, which points at the assumed void
    fraction rather than the correlation -- voidage is the one number here that
    was guessed rather than reported.
    """
    z, T, p, rho = _kesten_profile()
    mech = HydrazineShell405()
    G = 3.0 * LB / FT ** 2
    bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)

    i, j = -2, -1
    measured = p[i] - p[j]
    mu = mech.viscosity(0.5 * (T[i] + T[j]))
    predicted = 0.5 * (bed.pressure_gradient(z[i], G, rho[i], mu)
                       + bed.pressure_gradient(z[j], G, rho[j], mu)) * (z[j] - z[i])

    assert measured / PSI == pytest.approx(7.18, abs=0.1)      # guards the data
    assert 0.6 < predicted / measured < 1.4


@pytest.mark.skipif(not REFERENCE.exists(), reason="Kesten reference data absent")
def test_kesten_bed_is_inertia_dominated():
    """Particle Reynolds numbers are in the hundreds to thousands.

    Worth asserting because it justifies keeping Ergun's quadratic term: the
    gas accelerates enormously as it decomposes, so a Darcy-only model would
    be wrong by an order of magnitude by the bed exit.
    """
    z, T, p, rho = _kesten_profile()
    mech = HydrazineShell405()
    G = 3.0 * LB / FT ** 2
    bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)
    for zi, Ti, ri in zip(z, T, rho):
        m = bed.medium_at(zi)
        mu = mech.viscosity(Ti)
        assert m.reynolds(G, mu) > 100.0
        u = G / ri
        viscous = 150 * mu * (1 - m.eps) ** 2 / (m.eps ** 3 * m.d_p ** 2) * u
        inertial = 1.75 * ri * (1 - m.eps) / (m.eps ** 3 * m.d_p) * u * u
        assert inertial > 5.0 * viscous
