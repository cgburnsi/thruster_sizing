import math
import numpy as np

def _calc_diffusivity(T, P, ref_diff, ref_P=14.7, ref_T=492.0):
    """
    Replaces Fortran line 380: DI3 = D03*14.7/P*(T/492.)**1.823
    Adjusts reference diffusivity for Temperature and Pressure.
    """
    if P <= 0: return ref_diff
    # Power law scaling for gas diffusivity
    return ref_diff * (ref_P / P) * ((T / ref_T) ** 1.823)

def _calc_mass_transfer_coeff(G, rho, mu, diff, AP):
    """
    Replaces Fortran function KCF (Line 180).
    Calculates Mass Transfer Coefficient (Kc).
    Formula: .61 * (G/rho) * Sc^(-2/3) * Re^(-0.41)
    """
    if rho <= 0 or diff <= 0 or mu <= 0: return 0.0
    
    schmidt = mu / (rho * diff)
    reynolds = G / (AP * mu)
    
    # Avoid div by zero in fractional powers if Re is tiny
    if reynolds < 1e-9: reynolds = 1e-9
    
    # Fortran: .61*A/B*(C/(B*D))**-.667*(A/(E*C))**-.41
    # A=G, B=Rho, C=Mu, D=Diff, E=AP
    # Term 2: (Mu / (Rho*Diff))^-0.667 -> Sc^-0.667
    # Term 3: (G / (AP*Mu))^-0.41      -> Re^-0.41
    
    return 0.61 * (G / rho) * (schmidt ** -0.667) * (reynolds ** -0.41)

def _calc_heat_transfer_coeff(G, cp_bar, mu, AP):
    """
    Replaces Fortran line 470.
    Calculates Heat Transfer Coefficient (Hc).
    Formula: 0.74 * G * Cp * Re^(-0.41)
    """
    if mu <= 0 or AP <= 0: return 0.0
    reynolds = G / (AP * mu)
    if reynolds < 1e-9: reynolds = 1e-9
    
    # Fortran: .74*G*CFBAR*(G/(AP*MU))**-.41
    return 0.74 * G * cp_bar * (reynolds ** -0.41)

def _reaction_rate(conc, temp, k_arrhenius, n_order, beta, gamma, c_bulk):
    """
    Calculates local reaction rate incorporating the
    Langmuir-Hinshelwood / Arrhenius term from Fortran line 1670.
    """
    if conc <= 0: return 0.0
    
    # Normalized concentration
    ratio = conc / c_bulk
    if ratio < 0: ratio = 0
    
    # Denominator for heat generation term (1 + beta*(1-C/Cbulk))
    denom = 1.0 + beta * (1.0 - ratio)
    if abs(denom) < 1e-9: denom = 1e-9
    
    # Exponential Argument: Gamma * Beta * (1-C/Cbulk) / Denom
    exp_arg = gamma * beta * (1.0 - ratio) / denom
    
    # Rate = k * C^n * exp(...)
    # Note: The original code uses k0 * CI3^(1-n) * C^n
    # This effectively normalizes the rate constant units.
    pre_factor = k_arrhenius * (c_bulk ** (1.0 - n_order))
    
    return pre_factor * (conc ** n_order) * math.exp(exp_arg)

def _integrate_trapp(x0, x_end, n_steps, cps, ci3, tps, k, n, beta, gamma):
    """
    Replaces CALL TRAPP (Line 830).
    Integrates the reaction rate from penetration depth X0 to Surface (x=1.0).
    Assumes a linear concentration profile for the estimation.
    """
    if x_end <= x0: return 0.0
    
    dx = (x_end - x0) / n_steps
    integral = 0.0
    
    # Trapezoidal Rule
    # We iterate points. For linear profile:
    # Conc at X = CPS * (X - X0) / (1 - X0)
    
    x_points = np.linspace(x0, x_end, n_steps + 1)
    
    # Pre-calculate rates at all points
    rates = []
    for x in x_points:
        if x_end == x0:
            c_local = 0.0
        else:
            c_local = cps * (x - x0) / (x_end - x0)
            
        r = _reaction_rate(c_local, tps, k, n, beta, gamma, ci3)
        rates.append(r)
        
    # Summation (First + 2*Middle + Last) / 2 * dx
    integral = sum(rates) - 0.5 * (rates[0] + rates[-1])
    integral *= dx
    
    return integral

def SGRAD(T_bulk, G, concs_dict, A, AP, reaction, species, cfg, props, verbose=False, return_diag=False):
    """
    Solves the surface reaction/diffusion balance for a catalyst particle.
    
    Args:
        T_bulk: Bulk gas temperature [R]
        G: Bulk flow parameter [lb/ft^2-s]
        concs_dict: Dictionary of bulk concentrations [lb/ft^3]
        A: Catalyst geometric radius/characteristic length [ft]
        AP: External surface area per volume [1/ft]
        reaction: The Kinetics reaction object driving the gradient (e.g., R2/Ammonia)
        species: The Kinetics species object being consumed (e.g., NH3)
        cfg: ReactorConfig object
        props: PropertyTables object
        
    Returns:
        flux_mass: Mass flux of the target species [lb/ft^2-s]
        tgrad: Heat flux (technically h * deltaT) [BTU/ft^2-s]
        diag: Dictionary of internal solver states
    """
    
    # --- 1. Setup Parameters ---
    ci3 = concs_dict.get(species.name, 0.0)
    if ci3 <= 1e-12:
        if return_diag: return 0.0, 0.0, {}
        return 0.0, 0.0
        
    # Get species properties
    diff_ref = getattr(species, 'diffusivity', cfg.DIF3)
    
    # Get reaction kinetics
    # FIX: Access .A and .Ea_R directly (flattened), falling back to legacy config
    # Previously: k_arrhenius_A = reaction.kinetics.get("A", cfg.ALPHA2)
    k_arrhenius_A = getattr(reaction, "A", cfg.ALPHA2)
    
    # Previously: ea_R = reaction.kinetics.get("Ea_R", cfg.BGM)
    ea_R = getattr(reaction, "Ea_R", cfg.BGM)
    
    # Note: 'orders' is correctly accessed as a dictionary based on config.py
    n_order = reaction.orders.get(species.name, cfg.EN3)
    
    # Calculate Mixture Properties
    # Rho = Sum(Concentrations)
    rho_mix = sum(concs_dict.values())
    
    # Viscosity (Function of T)
    mu_mix = props.VISVST(T_bulk)
    
    # Specific Heats for calculating Cp_bar (Mixture)
    # CP_BAR = Sum(Ci * Cpi) / Rho
    cp_sum = 0.0
    for sp_name, conc in concs_dict.items():
        # Get Cp for this species. 
        # (Assuming props has a generic lookup or specific tables)
        # Mapping to legacy tables for simplicity based on names:
        val = 0.0
        if sp_name == "H2": val = props.CFTBL1(T_bulk)
        elif sp_name == "N2": val = props.CFTBL2(T_bulk)
        elif sp_name == "NH3": val = props.CFTBL3(T_bulk)
        elif sp_name == "N2H4": val = props.CFTBL4(T_bulk)
        cp_sum += conc * val
        
    cp_bar = cp_sum / rho_mix if rho_mix > 0 else 0.25
    
    # --- 2. Iteration Initialization (Fortran Lines 230-380) ---
    max_iter = 50
    tolerance = 0.05 # 5%
    
    # Initial guesses
    tps = T_bulk           # Surface Temp
    cps = ci3 / 2.0        # Surface Conc
    x0 = 0.0               # Penetration Depth
    
    # Damping factors
    waf1 = 0.8
    waf2 = 0.2
    
    # Loop variables
    converged = False
    flux_mass = 0.0
    flux_heat = 0.0
    
    # Diagnostics
    diag = {}

    for ii in range(max_iter):
        
        # A. Property Updates at Surface T (TPS)
        # Note: Fortran updates properties at T_bulk initially, then TPS.
        # We use TPS for the physics check.
        
        # Diffusivity (corrected for P and T)
        diff_eff = _calc_diffusivity(tps, cfg.PRES, diff_ref)
        
        # Mass Transfer Coeff (KC3)
        kc3 = _calc_mass_transfer_coeff(G, rho_mix, mu_mix, diff_eff, AP)
        
        # Heat Transfer Coeff (HC)
        # Note: Fortran uses average properties, here we simplify to local T
        hc = _calc_heat_transfer_coeff(G, cp_bar, mu_mix, AP)
        
        # Reaction Enthalpy (H3) at Surface T
        # (Assume Ammonia dissociation, lookup H3TBL)
        h_rxn = props.H3TBL(tps) 
        
        # Arrhenius Rate Constant at Surface T
        # k = A * exp(-Ea/T)
        k_rate = k_arrhenius_A * math.exp(-ea_R / tps)
        
        # Dimensionless Groups (Gamma, Beta)
        # Gamma = Ea / T
        gamma = ea_R / tps
        
        # Beta = - (Conc_Surface * Heat_Rxn * Diff) / (Cond * T) ???
        # Actually Fortran Line 660: BETA = -CPS*H3*DP3/(KP*TPS)
        # Wait, DP3 is Diffusivity term (14.7*D/P...).
        # KP is Thermal Conductivity? 
        # In Fortran COMMON/CO/ KP is 0.4e-4.
        
        # Let's use the explicit Fortran definition:
        # DP3 in Fortran is the Diffusivity term we calc'd as diff_eff
        beta_term = -cps * h_rxn * diff_eff / (cfg.KP * tps)
        
        # B. Calculate Gradient (DCPDX) from Mass Transfer Boundary
        # Flux = kc * (C_bulk - C_surface)
        # Gradient D(Conc)/D(x) = Flux / Diffusivity
        # In Fortran: DCPDX = KC3/DP3 * (CI3 - CPS)
        
        if diff_eff > 1e-12:
            dcpdx = (kc3 / diff_eff) * (ci3 - cps)
        else:
            dcpdx = 0.0
            
        # C. Calculate New Surface Temperature (Heat Balance)
        # Flux_Heat = HC * (T_bulk - T_surface)
        # Heat_Gen  = Flux_Mass * H_Rxn
        # Balance: HC*(T-TPS) = Flux_Mass * H_Rxn
        # TPS = T - (Flux_Mass * H_Rxn) / HC
        # Flux_Mass is approx (KC * (C - CPS))
        
        # Fortran Line 580: TPS = T - (H4*... + H3*DP3*DCPDX)/HC
        # Assuming only Ammonia (H3) is reacting here for SGRAD logic:
        flux_term = h_rxn * diff_eff * dcpdx
        
        if hc > 1e-9:
            tps_new = T_bulk - flux_term / hc
        else:
            tps_new = T_bulk
            
        # Sanity check T
        if tps_new < 100: tps_new = 100
        
        # D. Calculate Penetration Depth (X0) - Newton Raphson
        # Fortran Line 700: X0 = A - CPS/DCPDX
        # A is particle radius.
        
        if abs(dcpdx) > 1e-9:
            x0_calc = A - cps / dcpdx
        else:
            x0_calc = 0.0
            
        if x0_calc < 0: x0_calc = 0.0
        
        # Weighted average for X0 stability
        # Initialize x0 on first run
        if ii == 0: x0 = x0_calc
        x0 = waf1 * x0 + waf2 * x0_calc
        
        # E. Integrate Reaction Rate (Trapezoidal)
        riesum = _integrate_trapp(
            x0, A, 50, cps, ci3, tps, 
            k_rate, n_order, beta_term, gamma
        )
        
        # F. Update Surface Concentration (CPS)
        # Balance: Mass Transfer Flux = Total Reaction Rate inside particle
        # KC * (CI3 - CPS) = Integral_Rate
        # CPS = CI3 - Integral_Rate / KC
        
        if kc3 > 1e-12:
            cps_new = ci3 - riesum / kc3
        else:
            cps_new = ci3
            
        if cps_new < 0: cps_new = 0.0
        
        # G. Convergence Check
        # Check relative change in CPS and TPS
        diff_c = abs(cps - cps_new) / (ci3 + 1e-9)
        diff_t = abs(tps - tps_new) / (T_bulk + 1e-9)
        
        cps = cps_new
        tps = 0.5 * tps + 0.5 * tps_new # Damped T update
        
        if diff_c < tolerance and diff_t < tolerance:
            converged = True
            flux_mass = dcpdx * diff_eff # Gradient * Diffusivity
            flux_heat = hc * (T_bulk - tps)
            break
            
        # Adjust weights if struggling (Fortran Line 1390 logic)
        if ii > 25:
            waf1 = min(0.95, waf1 + 0.05)
            waf2 = 1.0 - waf1

    # Loop End
    if not converged and verbose:
        print(f"SGRAD Warning: Max iterations reached. Error C: {diff_c:.4f}, T: {diff_t:.4f}")
        
    # Final fluxes
    # Recalculate one last time to ensure consistency
    diff_eff = _calc_diffusivity(tps, cfg.PRES, diff_ref)
    if diff_eff > 0:
        dcpdx = (kc3 / diff_eff) * (ci3 - cps)
    flux_mass = dcpdx * diff_eff
    flux_heat = hc * (T_bulk - tps)
    
    if return_diag:
        diag = {
            "tps": tps,
            "cps": cps,
            "x0": x0,
            "kc": kc3,
            "hc": hc,
            "diff": diff_eff,
            "dcpdx": dcpdx,
            "iters": ii
        }
        return flux_mass, flux_heat, diag
        
    return flux_mass, flux_heat