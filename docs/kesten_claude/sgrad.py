import numpy as np

def SGRAD(
    TEMP,
    G,
    concentrations,
    A,
    AP,
    target_reaction,
    target_species,
    cfg,
    props,
    verbose=False,
    return_diag=False,
):
    """
    Compute diffusive and thermal fluxes at the catalyst surface.
    Generalized to solve for a specific target species and reaction.
    """

    # --- Small guards / constants ---
    eps = 1e-30
    NX = 24
    NPART = 50

    if A <= 0 or AP <= 0 or cfg.PRES <= 0 or TEMP <= 0:
        grad = 0.0
        tgrad = 0.0
        if return_diag:
            return grad, tgrad, {"note": "degenerate geometry / conditions"}
        return grad, tgrad

    # Bulk properties
    T = float(TEMP)
    P = float(cfg.PRES)

    # Calculate mixture density
    rho = sum(concentrations.values())
    if rho <= eps: rho = eps

    mu = float(props.VISVST(T))

    # Mixture Heat Capacity (cfbar)
    cfbar = 0.0
    for s_name, c_val in concentrations.items():
        sp = cfg.mechanism.get_species(s_name)
        if sp:
            cfbar += (c_val / rho) * sp.get_cp(T, props)

    # Diffusivity of the target species
    dp_bulk = _dp3f(T, target_species.diffusivity, P)

    # Mass transfer coefficient
    kc = _kcf(G, rho, mu, dp_bulk, AP)

    # Heat transfer coefficient
    hc = _hcf(G, mu, AP, cfbar)

    # Heat of Reaction for the target reaction
    h_rxn = target_species.get_enthalpy_reaction(T, props)
    
    # Target species concentration (bulk)
    ci_bulk = concentrations.get(target_species.name, 0.0)

    # Phase I (robust initial guess)
    x0a, cps, tps, dp_surf = _phase1_guess(
        T=T, P=P, A=A,
        ci_bulk=ci_bulk,
        kc=kc, hc=hc, h_rxn=h_rxn,
        target_reaction=target_reaction,
        target_species=target_species,
        concentrations=concentrations, 
        npart=NPART,
        props=props, verbose=verbose,
        cfg=cfg,
        # NEW: Pass fluid properties for correct k_c calculation
        G=G, rho=rho, mu=mu, AP=AP 
    )

    # Phase II (integral-equation iteration)
    diag = { "kc": kc, "hc": hc, "x0a": x0a }

    if ci_bulk <= 0.0:
        grad = 0.0
        tps = T 
        tgrad = hc * (T - tps)
        if return_diag:
            diag.update({"cps": 0.0, "tps": tps, "dp": dp_surf, "note": "ci<=0"})
            return grad, tgrad, diag
        return grad, tgrad

    x_nodes = np.linspace(x0a, 1.0, NX + 1)
    denom = max(1.0 - x0a, eps)
    cp_nodes = cps * (x_nodes - x0a) / denom
    cp_nodes[0] = 0.0

    relax_cp = 0.2
    prev_tmtp = None
    prev_cmcp = None

    for lp2 in range(1, 51):
        # Update properties at Surface Temperature (tps)
        h_rxn_s = target_species.get_enthalpy_reaction(tps, props)
        dp_surf = _dp3f(tps, target_species.diffusivity, P)

        gamma = target_reaction.Ea_R / max(tps, eps)
        
        beta = -cp_nodes[-1] * h_rxn_s * dp_surf / (float(cfg.KP) * max(tps, eps)) if cfg.KP != 0 else 0.0
        
        k0_base = target_reaction.A * np.exp(-gamma)
        
        k0 = k0_base
        for reac_sp, order in target_reaction.orders.items():
            if reac_sp != target_species.name:
                c_inhib = concentrations.get(reac_sp, 0.0)
                floor = getattr(cfg, "INHIBITOR_FLOOR", 1e-12)
                
                if order < 0.0:
                    c_eff = max(c_inhib, floor)
                    k0 *= (c_eff ** order)
                else:
                    k0 *= (max(c_inhib, 0.0) ** order)

        vnu = -kc / max(dp_surf, eps)
        avnu = A * vnu
        if abs(avnu) < 1e-12:
            ctrm = 1.0
        else:
            ctrm = (avnu + 1.0) / avnu

        cp_mid = 0.5 * (cp_nodes[:-1] + cp_nodes[1:])
        
        target_order = target_reaction.orders.get(target_species.name, 1.0)
        r_mid = _reaction_rate(cp_mid, ci_bulk, k0, gamma, beta, n_order=target_order)

        xl = x_nodes[:-1]
        xu = x_nodes[1:]
        e1 = _eval1(xl, xu)
        e2 = _eval2(xl, xu)

        dI1 = r_mid * e1
        dJ = r_mid * (e2 - ctrm * e1)

        I1 = np.concatenate(([0.0], np.cumsum(dI1)))
        I2 = np.zeros(NX + 1)
        I2[:-1] = np.cumsum(dJ[::-1])[::-1]

        fac = (A * A) / max(dp_surf, eps)
        cp_new = np.empty_like(cp_nodes)

        cp_new[0] = ci_bulk - fac * I2[0]
        for j in range(1, NX + 1):
            xj = x_nodes[j]
            cp_new[j] = ci_bulk - fac * ((1.0 / max(xj, eps) - ctrm) * I1[j] + I2[j])

        cp_new = np.clip(cp_new, 0.0, ci_bulk)

        cps_new = float(cp_new[-1])
        cmcp_new = ci_bulk - cps_new

        grad = kc * (ci_bulk - cps_new)

        # Estimate Tps (Simplified Thermal Balance for Iteration)
        term_other = 0.0
        for sp in cfg.mechanism.species:
            if sp.name != target_species.name and sp.name == "N2H4": 
                dp_other = _dp3f(tps, sp.diffusivity, P)
                kc_other = _kcf(G, rho, mu, dp_other, AP)
                h_other = sp.get_enthalpy_reaction(tps, props)
                c_other = concentrations.get(sp.name, 0.0)
                term_other += h_other * kc_other * c_other

        tps_guess = T - (term_other + h_rxn_s * grad) / max(hc, eps)
        if tps_guess < 0.0: tps_guess = 1.0
        
        tps_new = float(tps_guess)
        tmtp_new = T - tps_new

        if verbose and (lp2 == 1 or lp2 % 5 == 0):
            print(f"  [SGRAD] it={lp2:02d} x0a={x0a:.4f} cps={cps_new:.3e} TPS={tps_new:.2f}")

        if prev_tmtp is not None and prev_cmcp is not None:
            ok_t = abs(prev_tmtp - tmtp_new) / max(abs(tmtp_new), 1e-12) < 0.05
            ok_c = abs(prev_cmcp - cmcp_new) / max(abs(cmcp_new), 1e-12) < 0.05
            if ok_t and ok_c:
                cp_nodes = cp_new
                tps = tps_new
                break

        prev_tmtp = tmtp_new
        prev_cmcp = cmcp_new

        cp_nodes = (1.0 - relax_cp) * cp_nodes + relax_cp * cp_new
        cp_nodes[0] = 0.0
        cp_nodes[-1] = float(cp_nodes[-1])
        tps = 0.5 * tps + 0.5 * tps_new

    grad = kc * (ci_bulk - float(cp_nodes[-1]))
    tgrad = hc * (T - float(tps))

    if return_diag:
        diag.update({
            "cps": float(cp_nodes[-1]),
            "tps": float(tps),
            "h_rxn": float(target_species.get_enthalpy_reaction(tps, props)),
        })
        return grad, tgrad, diag

    return grad, tgrad


# ------------------------- helpers -------------------------

def _dp3f(T, D0, P):
    T, P, D0 = float(T), float(P), float(D0)
    if T <= 0 or P <= 0 or D0 <= 0: return 0.0
    base = (14.7 * D0 / P) * (T / 492.0) ** 1.823
    corr = 1.0 - np.exp(-0.0672 * P * 492.0 / (14.7 * T))
    return float(base * corr)

def _kcf(G, rho, mu, d, ap):
    G, rho, mu, d, ap = float(G), float(rho), float(mu), float(d), float(ap)
    if rho <= 0 or mu <= 0 or d <= 0 or ap <= 0: return 0.0
    return float(0.61 * (G / rho) * (mu / (rho * d)) ** (-0.667) * (G / (ap * mu)) ** (-0.41))

def _hcf(G, mu, ap, cfbar):
    G, mu, ap, cfbar = float(G), float(mu), float(ap), float(cfbar)
    if mu <= 0 or ap <= 0: return 0.0
    return float(0.74 * G * cfbar * (G / (ap * mu)) ** (-0.41))

def _reaction_rate(cp, ci_bulk, k0, gamma, beta, n_order=1.0):
    cp = np.asarray(cp, dtype=float)
    ci_bulk, k0, gamma, beta, n = float(ci_bulk), float(k0), float(gamma), float(beta), float(n_order)
    if ci_bulk <= 0 or k0 == 0: return np.zeros_like(cp)
    
    cp_pos = np.maximum(cp, 0.0)
    with np.errstate(over="ignore", invalid="ignore"):
        fact1 = k0 * (ci_bulk ** (1.0 - n)) * (cp_pos ** n)
        
    theta = 1.0 - (cp_pos / max(ci_bulk, 1e-30))
    denom = 1.0 + beta * theta
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    expo = np.clip(gamma * beta * theta / denom, -100.0, 100.0)
    return fact1 * np.exp(expo)

def _trapp_x2_rhet(x0a, cps, ci_bulk, k0, gamma, beta, n_order, npart):
    x0a, cps = float(x0a), float(cps)
    if x0a >= 1.0: return 0.0
    x = np.linspace(x0a, 1.0, int(npart) + 1)
    denom = max(1.0 - x0a, 1e-30)
    cp = cps * (x - x0a) / denom
    cp[0] = 0.0
    r = _reaction_rate(cp, ci_bulk, k0, gamma, beta, n_order=n_order)
    y = (x * x) * r
    return float(np.trapz(y, x))

def _phase1_guess(
    T, P, A, 
    ci_bulk, kc, hc, h_rxn, 
    target_reaction, target_species, 
    concentrations, npart, props, verbose, 
    cfg,
    G, rho, mu, AP # <--- ADDED ARGS HERE
):
    eps = 1e-30
    
    tps = T 
    cps = 0.5 * max(ci_bulk, 0.0)
    x0a = 0.0
    dp_surf = _dp3f(tps, target_species.diffusivity, P)

    for it in range(1, 31):
        h_rxn_s = target_species.get_enthalpy_reaction(tps, props)
        dp_surf = _dp3f(tps, target_species.diffusivity, P)
        
        grad = kc * (ci_bulk - cps)
        
        # Estimate Other Reaction Flux (N2H4) for Heat Balance
        term_other = 0.0
        for sp in cfg.mechanism.species:
            if sp.name != target_species.name and sp.name == "N2H4":
                dp_other = _dp3f(tps, sp.diffusivity, P)
                # Use real properties for kc calculation
                kc_other = _kcf(G, rho, mu, dp_other, AP)
                h_other = sp.get_enthalpy_reaction(tps, props)
                
                # Assume other concentration is roughly bulk (Phase 1 simplicity)
                c_other = concentrations.get(sp.name, 0.0)
                term_other += h_other * kc_other * c_other

        tps_new = T - (term_other + h_rxn_s * grad) / max(hc, eps)
        if tps_new < 0.0: tps_new = 1.0

        dcpdx = kc / max(dp_surf, eps) * (ci_bulk - cps)
        if dcpdx > 0.0:
            x0 = A - cps / max(dcpdx, eps)
        else:
            x0 = 0.0
        x0 = min(max(x0, 0.0), A)
        x0a = x0 / A

        gamma = target_reaction.Ea_R / max(tps, eps)
        
        beta = -cps * h_rxn_s * dp_surf / (float(cfg.KP) * max(tps, eps)) if cfg.KP != 0 else 0.0
        
        k0 = target_reaction.A * np.exp(-gamma)
        for reac_sp, order in target_reaction.orders.items():
            if reac_sp != target_species.name:
                c_inhib = concentrations.get(reac_sp, 0.0)
                if c_inhib > 0: k0 *= (c_inhib ** order)

        target_order = target_reaction.orders.get(target_species.name, 1.0)
        riesum = _trapp_x2_rhet(x0a, cps, ci_bulk, k0, gamma, beta, target_order, npart)
        
        cps_new = ci_bulk - A * riesum / max(kc, eps)
        cps_new = float(np.clip(cps_new, 0.0, ci_bulk))

        cps = 0.2 * cps_new + 0.8 * cps
        tps = 0.5 * tps_new + 0.5 * tps

        if abs(cps_new - cps) / max(ci_bulk, 1e-12) < 0.01 and abs(tps_new - tps) / max(tps_new, 1e-12) < 0.01:
            break

    return float(x0a), float(cps), float(tps), float(dp_surf)

def _eval1(a, b):
    return (b ** 3 - a ** 3) / 3.0

def _eval2(a, b):
    return (b ** 2 - a ** 2) / 2.0