"""Quasi-1-D isentropic reference solution.

Used for two things: initialising the 2-D field (which cuts the transient
dramatically and keeps startup robust), and providing the inviscid,
ideal-expansion benchmark that the converged CFD result is measured against.
"""
import numpy as np

from .thermo import G0


def mach_from_area_ratio(AR, gamma, supersonic=True, iters=100):
    """Invert the isentropic area-Mach relation by bisection.

    Bisection rather than Newton because the relation is flat near M = 1 and
    steep at large area ratios; robustness matters more than iteration count
    here, and the whole array converges in one vectorised sweep.
    """
    AR = np.asarray(AR, dtype=float)
    g = gamma
    ex = (g + 1.0) / (2.0 * (g - 1.0))

    def f(M):
        return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ex

    if supersonic:
        lo = np.full_like(AR, 1.0)
        hi = np.full_like(AR, 100.0)
    else:
        lo = np.full_like(AR, 1e-6)
        hi = np.full_like(AR, 1.0)

    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        val = f(mid)
        if supersonic:
            # f increases with M on the supersonic branch
            too_small = val < AR
            lo = np.where(too_small, mid, lo)
            hi = np.where(too_small, hi, mid)
        else:
            # f decreases with M on the subsonic branch
            too_big = val > AR
            lo = np.where(too_big, mid, lo)
            hi = np.where(too_big, hi, mid)
    return 0.5 * (lo + hi)


def isentropic_field(x, contour, gas, p0, T0):
    """Isentropic quasi-1-D state at each axial station in ``x``."""
    A = contour.area(x)
    At = np.pi * contour.r_throat ** 2
    AR = np.maximum(A / At, 1.0 + 1e-12)
    xt = contour.x_throat

    M = np.where(x <= xt,
                 mach_from_area_ratio(AR, gas.gamma, supersonic=False),
                 mach_from_area_ratio(AR, gas.gamma, supersonic=True))
    T = T0 / gas.T_ratio(M)
    p = p0 / gas.p_ratio(M)
    rho = p / (gas.R * T)
    a = np.sqrt(gas.gamma * gas.R * T)
    return dict(x=x, A=A, M=M, T=T, p=p, rho=rho, u=M * a, a=a)


def ideal_performance(contour, gas, p0, T0, p_amb):
    """Closed-form inviscid performance for the given contour.

    This is the same physics as ``thruster_sizing.py``, arranged as a forward
    problem (given p0, find thrust) instead of an inverse one.
    """
    At = np.pi * contour.r_throat ** 2
    Ae = np.pi * contour.r_exit ** 2
    eps = Ae / At
    Me = float(mach_from_area_ratio(np.array(eps), gas.gamma, supersonic=True))
    Te = T0 / gas.T_ratio(Me)
    pe = p0 / gas.p_ratio(Me)
    ue = Me * np.sqrt(gas.gamma * gas.R * Te)

    cstar = gas.cstar_ideal(T0)
    mdot = p0 * At / cstar
    F = mdot * ue + (pe - p_amb) * Ae
    return dict(M_exit=Me, T_exit=Te, p_exit=pe, u_exit=ue, area_ratio=eps,
                cstar=cstar, mdot=mdot, thrust=F, Isp=F / (mdot * G0),
                Cf=F / (p0 * At), throat_area=At, exit_area=Ae)


def initial_field(grid, gas, p0, T0, p_amb, ramp=0.0):
    """Build a starting primitive field from the quasi-1-D solution.

    ``ramp`` blends toward uniform stagnation conditions; 0 uses the full
    quasi-1-D solution, 1 starts the whole domain at rest.
    """
    sol = isentropic_field(grid.xc, grid.contour, gas, p0, T0)
    rho = sol["rho"]
    p = np.maximum(sol["p"], 0.5 * p_amb)
    u = sol["u"]
    v = np.zeros_like(u)

    if ramp > 0.0:
        rho0 = p0 / (gas.R * T0)
        rho = (1 - ramp) * rho + ramp * rho0
        p = (1 - ramp) * p + ramp * p0
        u = (1 - ramp) * u

    # Give the radial velocity the local wall slope so the initial field is
    # roughly aligned with the contour instead of purely axial.
    rw = grid.contour.r_wall(grid.xc)
    h = 1e-6 * (grid.contour.x_end - grid.contour.x_start)
    drdx = (grid.contour.r_wall(grid.xc + h) - grid.contour.r_wall(grid.xc - h)) / (2 * h)
    v = u * drdx * (grid.rc / np.maximum(rw, 1e-30))

    return np.stack([rho, u, v, p])
