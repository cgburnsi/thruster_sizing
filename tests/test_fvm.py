"""Verification tests for the axisymmetric finite-volume solver.

These are ordered from the smallest testable pieces (thermodynamics, grid
metrics, flux functions) up to whole-solver behaviour, so a failure points at
a layer rather than at "the CFD is wrong".
"""
import numpy as np
import pytest

from fvm import (PerfectGas, LOX_LH2_OF5, ConicalNozzle, BellNozzle, Grid,
                 BoundaryConditions, NozzleSolver, quasi1d, post)
from fvm.riemann import FLUX_SCHEMES, roe_flux, hllc_flux, ausm_plus_up_flux
from fvm.reconstruct import LIMITERS, muscl, get_limiter
from fvm.viscous import _green_gauss


@pytest.fixture(scope="module")
def gas():
    return PerfectGas(**LOX_LH2_OF5)


@pytest.fixture(scope="module")
def contour():
    return ConicalNozzle(4.15e-3, 0.145e-3, 100.0, 10e-3)


@pytest.fixture(scope="module")
def grid(contour):
    return Grid.from_contour(contour, ni=60, nj=24, wall_spacing=0.02)


# ---------------------------------------------------------------- thermo ---
def test_gas_constant_and_specific_heats(gas):
    assert gas.R == pytest.approx(8314.462618 / 11.8)
    assert gas.cp - gas.cv == pytest.approx(gas.R)
    assert gas.cp / gas.cv == pytest.approx(gas.gamma)


def test_area_ratio_unity_at_sonic(gas):
    assert gas.area_ratio(1.0) == pytest.approx(1.0)


def test_isentropic_relations_consistent(gas):
    M = np.array([0.1, 0.5, 1.0, 3.0, 5.4])
    # p0/p = (rho0/rho) * (T0/T)
    assert np.allclose(gas.p_ratio(M), gas.rho_ratio(M) * gas.T_ratio(M))


def test_viscosity_power_law_at_reference(gas):
    assert gas.viscosity(gas.T_mu_ref) == pytest.approx(gas.mu_ref)
    assert gas.conductivity(300.0) == pytest.approx(
        gas.cp * gas.viscosity(300.0) / gas.Pr)


# ------------------------------------------------------------------ grid ---
def test_volume_matches_analytic_integral(grid, contour):
    xs = np.linspace(contour.x_start, contour.x_end, 200001)
    exact = np.trapezoid(contour.r_wall(xs) ** 2 / 2.0, xs)
    assert grid.V.sum() == pytest.approx(exact, rel=2e-3)


def test_planar_area_matches_analytic_integral(grid, contour):
    xs = np.linspace(contour.x_start, contour.x_end, 200001)
    exact = np.trapezoid(contour.r_wall(xs), xs)
    assert grid.A_planar.sum() == pytest.approx(exact, rel=2e-3)


def test_closed_surface_identities(grid):
    """Sum of n*S over a cell must vanish axially and equal A_planar radially.

    The radial identity is exactly what makes the axisymmetric pressure source
    term consistent: uniform pressure produces no net radial acceleration only
    if sum(n_r S) == A_planar to round-off.
    """
    sx = (grid.nx_i[1:] * grid.S_i[1:] - grid.nx_i[:-1] * grid.S_i[:-1]
          + grid.nx_j[:, 1:] * grid.S_j[:, 1:] - grid.nx_j[:, :-1] * grid.S_j[:, :-1])
    sr = (grid.nr_i[1:] * grid.S_i[1:] - grid.nr_i[:-1] * grid.S_i[:-1]
          + grid.nr_j[:, 1:] * grid.S_j[:, 1:] - grid.nr_j[:, :-1] * grid.S_j[:, :-1])
    scale = grid.S_i.max()
    assert np.abs(sx).max() < 1e-12 * scale
    assert np.abs(sr - grid.A_planar).max() < 1e-12 * grid.A_planar.max()


def test_axis_faces_have_zero_area(grid):
    """Faces on r = 0 must carry no flux; this is what makes the axis BC free."""
    assert np.allclose(grid.S_j[:, 0], 0.0)


def test_green_gauss_exact_for_linear_field(grid):
    a, b = 3.0, 5.0
    gx, gr = _green_gauss(a * grid.xf_i + b * grid.rf_i,
                          a * grid.xf_j + b * grid.rf_j, grid)
    assert np.abs(gx - a).max() < 1e-9
    assert np.abs(gr - b).max() < 1e-9


def test_wall_clustering_is_toward_the_wall(contour):
    g = Grid.from_contour(contour, ni=20, nj=40, wall_spacing=0.004)
    eta = g.r_n[0] / g.r_n[0, -1]
    d = np.diff(eta)
    assert d[-1] == pytest.approx(0.004, rel=1e-3)
    assert d[-1] < d[0], "cells must be finer at the wall than on the axis"
    assert np.all(d > 0)


def test_bell_contour_builds_and_is_monotone_in_divergent():
    c = BellNozzle(4.15e-3, 0.145e-3, 100.0, 10e-3)
    x, r = c.polyline()
    div = x > c.x_arc_dn
    assert np.all(np.diff(r[div]) >= -1e-12)
    # exact because the sampler forces x = 0 (the true throat) into the polyline
    assert c.r_throat == pytest.approx(0.145e-3, rel=1e-12)
    assert c.area_ratio == pytest.approx(100.0, rel=1e-9)


# ------------------------------------------------------------- quasi-1-D ---
def test_mach_from_area_ratio_inverts_area_relation(gas):
    for branch in (True, False):
        M = np.array([0.05, 0.3, 0.8]) if not branch else np.array([1.5, 3.0, 5.4])
        AR = gas.area_ratio(M)
        assert np.allclose(quasi1d.mach_from_area_ratio(AR, gas.gamma, branch), M,
                           rtol=1e-6)


def test_ideal_performance_matches_hand_calculation(gas, contour):
    ideal = quasi1d.ideal_performance(contour, gas, 845e3, 3250.0, 666.6)
    assert ideal["M_exit"] == pytest.approx(5.4125, rel=1e-3)
    assert ideal["cstar"] == pytest.approx(gas.cstar_ideal(3250.0))
    assert ideal["mdot"] == pytest.approx(845e3 * ideal["throat_area"] / ideal["cstar"])


# ----------------------------------------------------------------- flux ----
def _uniform(rho=0.3, u=900.0, v=-120.0, p=4.0e5):
    """A single state shaped (4, 1, 1), matching the (4, ni, nj) face layout."""
    return np.array([rho, u, v, p]).reshape(4, 1, 1)


def _normal(nx, nr):
    return np.array([[nx]]), np.array([[nr]])


@pytest.mark.parametrize("name", sorted(FLUX_SCHEMES))
def test_flux_consistency(gas, name):
    """F(W, W, n) must reduce to the exact physical flux."""
    f = FLUX_SCHEMES[name]
    W = _uniform()
    nx, nr = _normal(0.6, 0.8)
    rho, u, v, p = W
    qn = u * nx + v * nr
    H = gas.gamma / (gas.gamma - 1) * p / rho + 0.5 * (u * u + v * v)
    exact = np.stack([rho * qn, rho * qn * u + p * nx,
                      rho * qn * v + p * nr, rho * qn * H])
    assert np.allclose(f(W, W, nx, nr, gas), exact, rtol=1e-10)


@pytest.mark.parametrize("name", sorted(FLUX_SCHEMES))
def test_flux_conservation_symmetry(gas, name):
    """Reversing the face must negate the flux: F(L,R,n) == -F(R,L,-n)."""
    f = FLUX_SCHEMES[name]
    WL = _uniform(0.30, 900.0, -120.0, 4.0e5)
    WR = _uniform(0.22, 1400.0, 60.0, 2.1e5)
    nx, nr = _normal(0.6, 0.8)
    a = f(WL, WR, nx, nr, gas)
    b = f(WR, WL, -nx, -nr, gas)
    assert np.allclose(a, -b, rtol=1e-10, atol=1e-8 * np.abs(a).max())


@pytest.mark.parametrize("name", sorted(FLUX_SCHEMES))
def test_flux_supersonic_is_full_upwind(gas, name):
    """With both states supersonic in +n, the flux must be the left flux."""
    f = FLUX_SCHEMES[name]
    WL = _uniform(0.10, 4000.0, 0.0, 2.0e4)
    WR = _uniform(0.08, 4300.0, 0.0, 1.2e4)
    nx, nr = _normal(1.0, 0.0)
    rho, u, v, p = WL
    H = gas.gamma / (gas.gamma - 1) * p / rho + 0.5 * u * u
    exact = np.stack([rho * u, rho * u * u + p, np.zeros_like(u), rho * u * H])
    assert np.allclose(f(WL, WR, nx, nr, gas), exact, rtol=1e-6)


# ------------------------------------------------------------- limiters ----
@pytest.mark.parametrize("name", sorted(LIMITERS))
def test_limiter_vanishes_at_extrema(name):
    """A sign change between neighbouring slopes must kill the reconstruction."""
    lim = get_limiter(name)
    a = np.array([1.0, -2.0, 3.0])
    b = np.array([-1.0, 3.0, -0.5])
    assert np.allclose(lim(a, b), 0.0, atol=1e-12)


@pytest.mark.parametrize("name", ["minmod", "vanleer", "vanalbada", "mc"])
def test_limiter_bounded_by_neighbour_slopes(name):
    lim = get_limiter(name)
    rng = np.random.default_rng(0)
    a = rng.uniform(0.1, 5.0, 500)
    b = rng.uniform(0.1, 5.0, 500)
    s = lim(a, b)
    assert np.all(s <= 2.0 * np.minimum(a, b) + 1e-12)
    assert np.all(s >= 0.0)


def test_muscl_reproduces_a_linear_field():
    """On a linear field the two face states must agree and hit the true value."""
    ng, n = 2, 12
    x = np.arange(-ng, n + ng, dtype=float)
    W = np.stack([1.0 + 0.1 * x, 2.0 - 0.3 * x, 0.5 * x, 3.0 + 0.2 * x])
    W = W[:, :, None]
    WL, WR = muscl(W, 1, ng, n, get_limiter("vanalbada"))
    assert np.allclose(WL, WR, rtol=1e-10)


# --------------------------------------------------------------- solver ----
def _solver(grid, gas, viscous, flux="roe", limiter="vanalbada"):
    bcs = BoundaryConditions(p0=845e3, T0=3250.0, p_amb=666.6)
    return NozzleSolver(grid, gas, bcs, flux=flux, limiter=limiter,
                        viscous=viscous, cfl=1.2, cfl_start=0.1, cfl_ramp=200)


def _uniform_residual(solver, state):
    W = np.empty_like(solver.W)
    for k, val in enumerate(state):
        W[k] = val
    return solver.residual_from_W(W)


def _cylindrical_grid(radius=2.0e-3, length=8.0e-3, ni=40, nj=20):
    x = np.linspace(0.0, length, ni + 1)
    eta = np.linspace(0.0, 1.0, nj + 1) ** 1.4      # deliberately non-uniform
    return Grid(np.repeat(x[:, None], nj + 1, axis=1),
                radius * np.repeat(eta[None, :], ni + 1, axis=0))


@pytest.mark.parametrize("state", [(0.30, 1500.0, 0.0, 4.0e5),
                                   (0.05, 3200.0, 0.0, 1.0e4)])
def test_free_stream_preservation_in_a_duct(gas, state):
    """Uniform axial flow in a constant-radius duct must give zero residual.

    This is the sharpest single check on the inviscid operator. It holds only
    if the metric identities close AND the axisymmetric pressure source is
    integrated over the planar cell area rather than the r-weighted volume --
    get that wrong and the radial momentum equation picks up a spurious
    acceleration proportional to p.
    """
    g = _cylindrical_grid()
    bcs = BoundaryConditions(p0=845e3, T0=3250.0, p_amb=666.6)
    s = NozzleSolver(g, gas, bcs, viscous=False)
    R = _uniform_residual(s, state)
    scale = max(state[0] * state[1] ** 2, state[3]) * g.S_i.max()
    assert np.abs(R).max() < 1e-12 * scale


def test_free_stream_preservation_on_curved_grid(grid, gas):
    """Same check on the real nozzle grid, away from the wall.

    The wall-adjacent layer is excluded on purpose: uniform axial flow *does*
    cross a sloped wall, and the solver correctly refuses to let it, so a
    non-zero residual there is the right answer rather than a metric error.
    """
    s = _solver(grid, gas, viscous=False)
    state = (0.30, 1500.0, 0.0, 4.0e5)
    R = _uniform_residual(s, state)[:, :, :-1]
    scale = max(state[0] * state[1] ** 2, state[3]) * grid.S_i.max()
    # ~1e-11 relative is round-off accumulated over the face sums; any genuine
    # metric or source-term error shows up at 1e-3 relative or worse.
    assert np.abs(R).max() < 1e-11 * scale


def test_conserved_primitive_roundtrip(grid, gas):
    s = _solver(grid, gas, viscous=True)
    rng = np.random.default_rng(1)
    W = np.stack([rng.uniform(0.01, 1.0, (grid.ni, grid.nj)),
                  rng.uniform(-500, 3000, (grid.ni, grid.nj)),
                  rng.uniform(-200, 200, (grid.ni, grid.nj)),
                  rng.uniform(1e3, 9e5, (grid.ni, grid.nj))])
    assert np.allclose(s.cons_to_prim(s.prim_to_cons(W)), W, rtol=1e-12)


def test_inviscid_wall_is_slip_and_viscous_wall_is_no_slip(grid, gas):
    from fvm import bc
    for viscous in (True, False):
        s = _solver(grid, gas, viscous=viscous)
        s.initialize(quasi1d.initial_field(grid, gas, 845e3, 3250.0, 666.6))
        W = s._fill(s.U)
        ng, nj = s.ng, grid.nj
        u_in, v_in = W[1, ng:-ng, ng + nj - 1], W[2, ng:-ng, ng + nj - 1]
        u_gh, v_gh = W[1, ng:-ng, ng + nj], W[2, ng:-ng, ng + nj]
        nx, nr = grid.nx_j[:, -1], grid.nr_j[:, -1]
        if viscous:
            assert np.allclose(0.5 * (u_in + u_gh), 0.0, atol=1e-9)
            assert np.allclose(0.5 * (v_in + v_gh), 0.0, atol=1e-9)
        else:
            # tangency: the face-averaged normal velocity vanishes, but the
            # tangential component is carried through unchanged
            qn = 0.5 * ((u_in + u_gh) * nx + (v_in + v_gh) * nr)
            assert np.allclose(qn, 0.0, atol=1e-6)
            tin = -u_in * nr + v_in * nx
            tgh = -u_gh * nr + v_gh * nx
            assert np.allclose(tin, tgh, rtol=1e-9)


def test_axis_ghost_cells_mirror_radial_velocity(grid, gas):
    s = _solver(grid, gas, viscous=True)
    s.initialize(quasi1d.initial_field(grid, gas, 845e3, 3250.0, 666.6))
    W = s._fill(s.U)
    ng = s.ng
    assert np.allclose(W[2, :, ng - 1], -W[2, :, ng])
    assert np.allclose(W[0, :, ng - 1], W[0, :, ng])
    assert np.allclose(W[3, :, ng - 1], W[3, :, ng])


@pytest.mark.slow
def test_inviscid_nozzle_approaches_quasi_1d(gas, contour):
    """An Euler run must recover the isentropic answer to a few per cent."""
    g = Grid.from_contour(contour, ni=140, nj=44, wall_spacing=0.02)
    s = _solver(g, gas, viscous=False, limiter="minmod")
    s.initialize(quasi1d.initial_field(g, gas, 845e3, 3250.0, 666.6))
    s.run(max_iter=5000, tol=1e-9, print_every=0)

    ideal = quasi1d.ideal_performance(contour, gas, 845e3, 3250.0, 666.6)
    st = post.station_integrals(s)
    # Mass flow through the choked throat is the least grid-sensitive quantity
    assert st["mdot"][-1] / ideal["mdot"] == pytest.approx(1.0, abs=0.03)
    assert s.fields()["M"][-1, 0] == pytest.approx(ideal["M_exit"], rel=0.05)

    # Mass must be conserved through the supersonic divergent section. The
    # window excludes the throat face itself (a face-averaged state is a poor
    # quadrature where the flow is transonic) and the last two faces (the
    # outflow stencil). The chamber is excluded too: its convective time scale
    # is ~800x the nozzle's, so it is still drifting at this iteration count.
    Ld = contour.x_end - contour.x_throat
    win = ((st["x"] > contour.x_throat + 0.15 * Ld)
           & (st["x"] < contour.x_throat + 0.90 * Ld))
    md = st["mdot"][win]
    assert md.size > 10
    assert (md.max() - md.min()) / md.mean() < 0.01


@pytest.mark.slow
def test_viscous_run_loses_thrust_and_stays_bounded(gas, contour):
    """The physical signature of a micro-nozzle: large, one-signed viscous loss."""
    g = Grid.from_contour(contour, ni=120, nj=48, wall_spacing=0.006)
    s = _solver(g, gas, viscous=True)
    s.initialize(quasi1d.initial_field(g, gas, 845e3, 3250.0, 666.6))
    s.run(max_iter=4000, tol=1e-9, print_every=0)

    ideal = quasi1d.ideal_performance(contour, gas, 845e3, 3250.0, 666.6)
    perf = post.performance(s, ideal)
    assert np.isfinite(s.U).all()
    assert 0.5 < perf["eta_thrust"] < 1.0, "viscous thrust must fall below ideal"
    assert 0.7 < perf["Cd"] < 1.0, "discharge coefficient must fall below unity"
    assert post.throat_reynolds(s) == pytest.approx(1170.0, rel=0.25)
    assert perf["y_plus_max"] < 5.0, "wall grid too coarse to resolve the layer"
