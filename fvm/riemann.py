"""Approximate Riemann solvers for the 2-D compressible Euler fluxes.

Every solver takes left/right primitive states ``W = [rho, u, v, p]`` as
``(4, ...)`` arrays together with the face unit normal ``(nx, nr)``, and
returns the normal flux ``(4, ...)`` of ``[rho, rho*u, rho*v, rho*E]``.
The arrays are shape-agnostic, so the same code serves both face families.
"""
import numpy as np

TINY = 1e-30


def _physical_flux(rho, u, v, p, nx, nr, gamma):
    qn = u * nx + v * nr
    H = gamma / (gamma - 1.0) * p / rho + 0.5 * (u * u + v * v)
    m = rho * qn
    return np.stack([m, m * u + p * nx, m * v + p * nr, m * H]), qn, H


def roe_flux(WL, WR, nx, nr, gas, entropy_fix=0.1):
    """Roe's flux-difference splitting with a Harten-Hyman entropy fix.

    The entropy fix is applied to the acoustic waves only; smearing the
    entropy and shear waves would thicken the boundary layer.
    """
    g = gas.gamma
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR

    FL, qnL, HL = _physical_flux(rL, uL, vL, pL, nx, nr, g)
    FR, qnR, HR = _physical_flux(rR, uR, vR, pR, nx, nr, g)

    sq = np.sqrt(rR / rL)
    d = 1.0 / (1.0 + sq)
    u_ = (uL + sq * uR) * d
    v_ = (vL + sq * vR) * d
    H_ = (HL + sq * HR) * d
    r_ = np.sqrt(rL * rR)
    q2 = u_ * u_ + v_ * v_
    a2 = np.maximum((g - 1.0) * (H_ - 0.5 * q2), TINY)
    a_ = np.sqrt(a2)
    qn_ = u_ * nx + v_ * nr

    drho = rR - rL
    dp = pR - pL
    dqn = qnR - qnL
    du = uR - uL
    dv = vR - vL

    l1 = qn_ - a_
    l2 = qn_
    l4 = qn_ + a_
    delta = entropy_fix * (np.abs(qn_) + a_) + TINY
    al1 = np.where(np.abs(l1) < delta, (l1 * l1 + delta * delta) / (2.0 * delta), np.abs(l1))
    al4 = np.where(np.abs(l4) < delta, (l4 * l4 + delta * delta) / (2.0 * delta), np.abs(l4))
    al2 = np.abs(l2)

    w1 = al1 * (dp - r_ * a_ * dqn) / (2.0 * a2)     # acoustic, qn - a
    w2 = al2 * (drho - dp / a2)                       # entropy
    w4 = al4 * (dp + r_ * a_ * dqn) / (2.0 * a2)     # acoustic, qn + a
    ws = al2 * r_                                     # shear
    dut = du - dqn * nx
    dvt = dv - dqn * nr

    diss = np.stack([
        w1 + w2 + w4,
        w1 * (u_ - a_ * nx) + w2 * u_ + w4 * (u_ + a_ * nx) + ws * dut,
        w1 * (v_ - a_ * nr) + w2 * v_ + w4 * (v_ + a_ * nr) + ws * dvt,
        w1 * (H_ - a_ * qn_) + w2 * 0.5 * q2 + w4 * (H_ + a_ * qn_)
        + ws * (u_ * dut + v_ * dvt),
    ])
    return 0.5 * (FL + FR) - 0.5 * diss


def hllc_flux(WL, WR, nx, nr, gas, **_):
    """HLLC with Einfeldt (Roe-averaged) wave-speed estimates.

    Positivity preserving and free of the carbuncle instability, at the cost
    of slightly more contact smearing than Roe.
    """
    g = gas.gamma
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR

    FL, qnL, HL = _physical_flux(rL, uL, vL, pL, nx, nr, g)
    FR, qnR, HR = _physical_flux(rR, uR, vR, pR, nx, nr, g)
    aL = np.sqrt(g * pL / rL)
    aR = np.sqrt(g * pR / rR)
    EL = pL / ((g - 1.0) * rL) + 0.5 * (uL * uL + vL * vL)
    ER = pR / ((g - 1.0) * rR) + 0.5 * (uR * uR + vR * vR)

    sq = np.sqrt(rR / rL)
    d = 1.0 / (1.0 + sq)
    u_ = (uL + sq * uR) * d
    v_ = (vL + sq * vR) * d
    H_ = (HL + sq * HR) * d
    a_ = np.sqrt(np.maximum((g - 1.0) * (H_ - 0.5 * (u_ * u_ + v_ * v_)), TINY))
    qn_ = u_ * nx + v_ * nr

    SL = np.minimum(qnL - aL, qn_ - a_)
    SR = np.maximum(qnR + aR, qn_ + a_)

    dL = rL * (SL - qnL)
    dR = rR * (SR - qnR)
    SM = (pR - pL + dL * qnL - dR * qnR) / np.where(np.abs(dL - dR) < TINY, TINY, dL - dR)

    def star(rho, u, v, E, qn, p, S, F):
        fac = rho * (S - qn) / np.where(np.abs(S - SM) < TINY, TINY, S - SM)
        Us = np.stack([
            fac,
            fac * (u + (SM - qn) * nx),
            fac * (v + (SM - qn) * nr),
            fac * (E + (SM - qn) * (SM + p / (rho * np.where(np.abs(S - qn) < TINY, TINY, S - qn)))),
        ])
        U = np.stack([rho, rho * u, rho * v, rho * E])
        return F + S * (Us - U)

    FsL = star(rL, uL, vL, EL, qnL, pL, SL, FL)
    FsR = star(rR, uR, vR, ER, qnR, pR, SR, FR)

    F = np.where(SL >= 0.0, FL,
                 np.where(SM >= 0.0, FsL,
                          np.where(SR >= 0.0, FsR, FR)))
    return F


def ausm_plus_up_flux(WL, WR, nx, nr, gas, M_inf=0.2, Kp=0.25, Ku=0.75, sigma=1.0, **_):
    """Liou's AUSM+-up.

    Carries a low-Mach fix, which in principle suits a nozzle whose chamber
    sits at M ~ 1e-3 while its exit runs at M ~ 5. In practice it is the
    weakest of the three schemes for this geometry, and it is not the default.

    The ``M_inf`` cutoff must not be set near the chamber Mach number: the
    upwind switch on ``sign(M_face)`` flips the donor density back and forth
    in a near-stagnant chamber, and with ``M_inf = 0.05`` this grows without
    bound in about 170 iterations. Values of 0.2 and above are stable but
    dissipate enough to cost roughly 2% in predicted mass flow relative to
    Roe on the same grid.
    """
    g = gas.gamma
    rL, uL, vL, pL = WL
    rR, uR, vR, pR = WR

    qnL = uL * nx + vL * nr
    qnR = uR * nx + vR * nr
    HL = g / (g - 1.0) * pL / rL + 0.5 * (uL * uL + vL * vL)
    HR = g / (g - 1.0) * pR / rR + 0.5 * (uR * uR + vR * vR)

    a = 0.5 * (np.sqrt(g * pL / rL) + np.sqrt(g * pR / rR))
    ML = qnL / a
    MR = qnR / a

    Mbar2 = 0.5 * (qnL * qnL + qnR * qnR) / (a * a)
    M02 = np.minimum(1.0, np.maximum(Mbar2, M_inf * M_inf))
    M0 = np.sqrt(M02)
    fa = M0 * (2.0 - M0)

    def M1p(M):
        return 0.5 * (M + np.abs(M))

    def M1m(M):
        return 0.5 * (M - np.abs(M))

    def M2p(M):
        return 0.25 * (M + 1.0) ** 2

    def M2m(M):
        return -0.25 * (M - 1.0) ** 2

    beta = 1.0 / 8.0
    sup = np.abs(ML) >= 1.0
    M4p = np.where(sup, M1p(ML), M2p(ML) * (1.0 - 16.0 * beta * M2m(ML)))
    sup = np.abs(MR) >= 1.0
    M4m = np.where(sup, M1m(MR), M2m(MR) * (1.0 + 16.0 * beta * M2p(MR)))

    rho_avg = 0.5 * (rL + rR)
    Mp = -Kp / fa * np.maximum(1.0 - sigma * Mbar2, 0.0) * (pR - pL) / (rho_avg * a * a)
    Mf = M4p + M4m + Mp

    alpha = 3.0 / 16.0 * (-4.0 + 5.0 * fa * fa)
    sup = np.abs(ML) >= 1.0
    P5p = np.where(sup, M1p(ML) / np.where(np.abs(ML) < TINY, TINY, ML),
                   M2p(ML) * ((2.0 - ML) - 16.0 * alpha * ML * M2m(ML)))
    sup = np.abs(MR) >= 1.0
    P5m = np.where(sup, M1m(MR) / np.where(np.abs(MR) < TINY, TINY, MR),
                   M2m(MR) * ((-2.0 - MR) + 16.0 * alpha * MR * M2p(MR)))

    pu = -Ku * P5p * P5m * (rL + rR) * fa * a * (qnR - qnL)
    pf = P5p * pL + P5m * pR + pu

    mdot = a * Mf * np.where(Mf > 0.0, rL, rR)
    psi1 = np.ones_like(mdot)
    psi2 = np.where(Mf > 0.0, uL, uR)
    psi3 = np.where(Mf > 0.0, vL, vR)
    psi4 = np.where(Mf > 0.0, HL, HR)

    return np.stack([
        mdot * psi1,
        mdot * psi2 + pf * nx,
        mdot * psi3 + pf * nr,
        mdot * psi4,
    ])


FLUX_SCHEMES = {
    "roe": roe_flux,
    "hllc": hllc_flux,
    "ausm": ausm_plus_up_flux,
}


def get_flux(name):
    try:
        return FLUX_SCHEMES[name.lower()]
    except KeyError:
        raise ValueError(f"unknown flux scheme {name!r}; "
                         f"choose from {sorted(FLUX_SCHEMES)}") from None
