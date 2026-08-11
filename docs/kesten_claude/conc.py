def CONC(T, P, R, H, HF, mechanism, props, seed_conversion=1e-3, h2_mole_floor=1e-12):
    """
    Vapor-region inlet concentrations using the legacy Kesten relation.

    The original model initializes vapor composition from the vapor-entry
    enthalpy and the hydrazine latent term H4(T):

        XV = -(H - HF) / H4
        C4 = (P*WM4)/(R*T) * (1 - XV)/(1 + XV)
        C3 = (P*WM3)/(R*T) * XV/(1 + XV)
        C2 = (P*WM2)/(2*R*T) * XV/(1 + XV)
        C1 = (P*WM1)/(2*R*T) * XV/(1 + XV)

    Returns mass concentrations in mechanism.species order.
    """
    if T <= 0.0 or R == 0.0:
        return [0.0 for _ in mechanism.species]

    sp_h2 = mechanism.get_species("H2")
    sp_n2 = mechanism.get_species("N2")
    sp_nh3 = mechanism.get_species("NH3")
    sp_n2h4 = mechanism.get_species("N2H4")
    if not all((sp_h2, sp_n2, sp_nh3, sp_n2h4)):
        return [0.0 for _ in mechanism.species]

    h4 = props.H4TBL(T)
    if h4 == 0.0:
        return [0.0 for _ in mechanism.species]

    xv = -(H - HF) / h4
    denom = 1.0 + xv
    if abs(denom) < 1e-30:
        return [0.0 for _ in mechanism.species]

    c4 = ((P * sp_n2h4.mw) / (R * T)) * ((1.0 - xv) / denom)
    c3 = ((P * sp_nh3.mw) / (R * T)) * (xv / denom)
    c2 = ((P * sp_n2.mw) / (2.0 * R * T)) * (xv / denom)
    c1 = ((P * sp_h2.mw) / (2.0 * R * T)) * (xv / denom)

    lookup = {
        "H2": max(0.0, c1),
        "N2": max(0.0, c2),
        "NH3": max(0.0, c3),
        "N2H4": max(0.0, c4),
    }
    return [lookup.get(sp.name, 0.0) for sp in mechanism.species]
