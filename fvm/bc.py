"""Boundary conditions, applied by filling two layers of ghost cells.

Boundaries
----------
i = 0      subsonic inlet on total conditions (p0, T0), flow along the face
           normal, with the outgoing Riemann invariant extrapolated
i = ni     outflow: supersonic extrapolation, switching pointwise to a
           characteristic back-pressure condition wherever the local normal
           Mach number drops below one (which it does inside the exit-plane
           boundary layer)
j = 0      symmetry axis
j = nj     solid wall: no slip, adiabatic or isothermal
"""
import numpy as np

TINY = 1e-30


class BoundaryConditions:
    """Container for the boundary data.

    Parameters
    ----------
    p0, T0 : float
        Inlet total pressure [Pa] and total temperature [K].
    p_amb : float
        Ambient (back) pressure [Pa].
    wall : {'adiabatic'} or float
        Wall thermal condition; a float is taken as a fixed wall temperature [K].
    inlet : {'stagnation', 'pressure', 'riemann'}
        How the subsonic inlet closes the system. All three impose p0 and T0;
        they differ in which piece of information they take from the interior.

        ``stagnation`` (default) extrapolates the area-averaged normal
        velocity and imposes it uniformly. This is the only one of the three
        that survives a high contraction ratio, and it is what you want here.

        ``pressure`` extrapolates static pressure and reads the Mach number
        off the isentropic p0/p relation. Because ``u ~ sqrt(p0 - p)``, it
        amplifies pressure noise: at CR = 819 the chamber dynamic pressure is
        ~0.3 Pa against a 845 kPa static level, so sub-Pa cell-to-cell
        variation produces order-100% scatter in inlet velocity.

        ``riemann`` extrapolates the outgoing ``qn - a`` invariant. It is the
        conventional textbook choice but is numerically hopeless here: it
        recovers a ~1 m/s chamber velocity by differencing two quantities of
        order 13000 m/s, which leaves no significant digits. Use it only for
        chamber Mach numbers above roughly 0.05.
    """

    def __init__(self, p0, T0, p_amb, wall="adiabatic", inlet="stagnation"):
        self.p0 = float(p0)
        self.T0 = float(T0)
        self.p_amb = float(p_amb)
        self.wall = wall
        if inlet not in ("stagnation", "pressure", "riemann"):
            raise ValueError("inlet must be 'stagnation', 'pressure' or 'riemann'")
        self.inlet = inlet

    @property
    def wall_is_adiabatic(self):
        return isinstance(self.wall, str)

    @property
    def T_wall(self):
        return None if self.wall_is_adiabatic else float(self.wall)


def apply(W, grid, gas, bcs, ng=2, viscous=True):
    """Fill the ghost cells of ``W`` (4, ni+2ng, nj+2ng) in place.

    ``viscous`` selects the wall treatment: no slip when True, and a slip
    (tangency) wall when False. An Euler run given no-slip ghosts would
    manufacture a shear layer at the wall through the MUSCL stencil of the
    last interior face, even though the wall flux itself is imposed.
    """
    ni, nj = grid.ni, grid.nj
    g = gas.gamma
    R = gas.R
    js = slice(ng, ng + nj)

    # ---- inlet (i = 0) : subsonic, total conditions ---------------------
    nx = grid.nx_i[0, :]
    nr = grid.nr_i[0, :]
    rho_i, u_i, v_i, p_i = W[:, ng, js]
    a_i = np.sqrt(g * p_i / rho_i)
    qn_i = u_i * nx + v_i * nr

    if bcs.inlet == "stagnation":
        # Take the area-averaged normal velocity from the interior and apply
        # it uniformly. Extrapolating velocity is well conditioned (it is an
        # O(1 m/s) quantity carried across at face value), and averaging over
        # the plane stops cell-to-cell scatter from seeding a spurious jet.
        S = grid.S_i[0, :]
        q_ext = float(np.sum(qn_i * S) / max(np.sum(S), TINY))
        q_ext = max(q_ext, 0.0)
        T_est = max(bcs.T0 - 0.5 * q_ext * q_ext / gas.cp, 0.1 * bcs.T0)
        M_b = np.full_like(qn_i, min(q_ext / np.sqrt(g * R * T_est), 0.99))
    elif bcs.inlet == "pressure":
        # Extrapolate static pressure, then read M off p0/p. expm1 keeps the
        # significant digits of a near-unity pressure ratio.
        p_b = np.minimum(p_i, bcs.p0)
        x = np.log(bcs.p0 / p_b) * (g - 1.0) / g
        M_b = np.sqrt(np.maximum(2.0 / (g - 1.0) * np.expm1(x), 0.0))
        M_b = np.minimum(M_b, 0.99)
    else:
        Jm = qn_i - 2.0 * a_i / (g - 1.0)      # outgoing invariant (qn - a wave)
        a0sq = g * R * bcs.T0
        c = 2.0 / (g - 1.0)
        disc = Jm * Jm - (1.0 + c) * (0.5 * (g - 1.0) * Jm * Jm - a0sq)
        a_b = (-Jm + np.sqrt(np.maximum(disc, 0.0))) / (1.0 + c)
        q_b = np.maximum(Jm + c * a_b, 0.0)
        M_b = np.minimum(q_b / np.maximum(a_b, TINY), 0.99)

    T_b = bcs.T0 / gas.T_ratio(M_b)
    p_b = bcs.p0 / gas.p_ratio(M_b)
    rho_b = p_b / (R * T_b)
    q_b = M_b * np.sqrt(g * R * T_b)

    for m in range(ng):
        W[0, m, js] = rho_b
        W[1, m, js] = q_b * nx
        W[2, m, js] = q_b * nr
        W[3, m, js] = p_b

    # ---- outlet (i = ni) ------------------------------------------------
    nx = grid.nx_i[-1, :]
    nr = grid.nr_i[-1, :]
    rho_i, u_i, v_i, p_i = W[:, ng + ni - 1, js]
    a_i = np.sqrt(g * p_i / rho_i)
    qn_i = u_i * nx + v_i * nr
    supersonic = qn_i >= a_i

    p_o = bcs.p_amb
    rho_o = np.maximum(rho_i + (p_o - p_i) / (a_i * a_i), 1e-12 * rho_i)
    fac = (p_i - p_o) / (rho_i * a_i)
    u_o = u_i + nx * fac
    v_o = v_i + nr * fac

    rho_g = np.where(supersonic, rho_i, rho_o)
    u_g = np.where(supersonic, u_i, u_o)
    v_g = np.where(supersonic, v_i, v_o)
    p_g = np.where(supersonic, p_i, p_o)

    for m in range(ng):
        W[0, ng + ni + m, js] = rho_g
        W[1, ng + ni + m, js] = u_g
        W[2, ng + ni + m, js] = v_g
        W[3, ng + ni + m, js] = p_g

    # ---- axis (j = 0) : mirror about r = 0 ------------------------------
    for m in range(ng):
        src = ng + m
        dst = ng - 1 - m
        W[0, :, dst] = W[0, :, src]
        W[1, :, dst] = W[1, :, src]
        W[2, :, dst] = -W[2, :, src]
        W[3, :, dst] = W[3, :, src]

    # ---- wall (j = nj) : no slip (viscous) or tangency (inviscid) --------
    # The wall normal is only defined over the interior i range, so edge-pad
    # it to cover the i-ghost columns and keep the corners consistent.
    nxw = np.concatenate([np.full(ng, grid.nx_j[0, -1]), grid.nx_j[:, -1],
                          np.full(ng, grid.nx_j[-1, -1])])
    nrw = np.concatenate([np.full(ng, grid.nr_j[0, -1]), grid.nr_j[:, -1],
                          np.full(ng, grid.nr_j[-1, -1])])

    for m in range(ng):
        src = ng + nj - 1 - m
        dst = ng + nj + m
        p_w = W[3, :, src]
        W[3, :, dst] = p_w
        if viscous:
            W[1, :, dst] = -W[1, :, src]
            W[2, :, dst] = -W[2, :, src]
        else:
            qn = W[1, :, src] * nxw + W[2, :, src] * nrw
            W[1, :, dst] = W[1, :, src] - 2.0 * qn * nxw
            W[2, :, dst] = W[2, :, src] - 2.0 * qn * nrw
        if bcs.wall_is_adiabatic or not viscous:
            W[0, :, dst] = W[0, :, src]
        else:
            T_in = W[3, :, src] / (R * W[0, :, src])
            T_g = np.maximum(2.0 * bcs.T_wall - T_in, 0.2 * bcs.T_wall)
            W[0, :, dst] = p_w / (R * T_g)
    return W
