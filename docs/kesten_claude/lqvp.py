from param import PARAM
from print_data import print_lqvp_data
from state import ReactorState
from config import ReactorConfig
from properties import PropertyTables

# Kesten liquid-vapor reference data (as provided)
# Z values [ft]
z_list = [
    0.0007376, 0.0007383, 0.0007390, 0.0007396, 0.0007403, 0.0007410,
    0.0007417, 0.0007424, 0.0007431, 0.0007437, 0.0007444, 0.0007449,
]

# TEMP values [degR]
temp_list = [
    820.0, 820.0, 820.0, 820.0, 820.0, 820.0,
    820.0, 820.0, 820.0, 820.0, 820.0, 820.0,
]

# H values [BTU/lb]
h_list = [
    212.63, 259.63, 306.63, 353.63, 400.64, 447.64,
    494.64, 541.64, 588.65, 635.65, 682.65, 715.48,
]

# WFV values [-]
wfv_list = [
    0.0, 0.0934717, 0.1869430, 0.2804150, 0.3738870, 0.4673580,
    0.5608300, 0.6543020, 0.7477730, 0.8412450, 0.9347170, 1.0,
]


def _print_lv_row(j_global, j_ref, z, temp, h, wfv):
    """Print one row like the liquid section: index, computed line, reference line."""
    print(j_global)
    print_lqvp_data(z, temp, h, wfv)
    if 0 <= j_ref < len(z_list):
        print_lqvp_data(z_list[j_ref], temp_list[j_ref], h_list[j_ref], wfv_list[j_ref])


def LQVP(state: ReactorState, cfg: ReactorConfig, props: PropertyTables):
    """
    Simulate the Liquid-Vapor transition region.
    
    Now operates directly on the ReactorState arrays, populating Z, H, TEMP, WFV, etc.
    """
    
    # Local index into the Kesten LV reference lists (0..len-1)
    j_ref = 0
    
    # Start at the current index from state
    jj = state.ii

    # Entry conditions
    H = state.H[jj]
    Z_prev = state.Z[jj]  # Previous Z (or current Z at start)
    
    # We need the previous deriv/dhdz or current? 
    # In the original, DERIV[JJ] was initialized from arguments.
    # The liquid loop left state.DERIV[jj] populated, so we are good.
    
    # Initialize first point in LV region
    state.TEMP[jj] = cfg.TVAP
    
    denom = (cfg.HV - cfg.HL)
    wfv = (H - cfg.HL) / denom if denom != 0.0 else 0.0
    state.WFV[jj] = wfv
    
    if cfg.print_rows:
        _print_lv_row(jj, j_ref, state.Z[jj], cfg.TVAP, H, wfv)

    jj += 1
    j_ref += 1

    while jj < cfg.max_steps:
        # In LV region, temperature is held at TVAP
        state.TEMP[jj] = cfg.TVAP
        
        # Derivative is held constant (simplified physics in LV region)
        state.DERIV[jj] = state.DERIV[jj - 1]

        # Use previous Z for lookups
        Z_prev = state.Z[jj - 1]
        
        H4 = props.H4TBL(cfg.TVAP)
        AP = props.ZTBLAP(Z_prev)
        # A  = props.ZTBLA(Z_prev) # Not used in eqn below, but available

        # Calculate G at current Z
        if (Z_prev - cfg.Z0) < 0:
            G = cfg.G0 + cfg.FC * Z_prev
        else:
            G = cfg.GATZ0
        
        # Recalculate parameters (DPA needed)
        # We can pass 0.0 for CN2H4 since it's not primary driver here? 
        # Original code called PARAM with 1, 0, 0, 1... let's replicate exactly.
        # Original: PARAM(TVAP, Z_prev, 1, 0.0, 0.0, 1, ...)
        _, _, _, _, DPA = PARAM(
            cfg.TVAP, Z_prev, 1, 0.0, 0.0, 1,
            cfg.Z0, cfg.G0, cfg.FC, cfg.GATZ0,
            cfg.AGM, cfg.BGM, cfg.ALPHA1, cfg.ALPHA2, 
            cfg.DIF3, cfg.DIF4,
            cfg.PRES, cfg.KP, cfg.C1
        )

        # Compute DHDZ
        # DHDZ = -(H4 * DPA * AP * DERIV + FC * (H - HF)) / G
        # Note: using 'H' from previous step in the FC term
        val_dhdz = -(H4 * DPA * AP * state.DERIV[jj] + cfg.FC * (H - cfg.HF)) / G
        state.DHDZ[jj] = val_dhdz

        # Compute step size DZ
        # DZ = -H4 / (ENMX2 * DHDZ)
        if val_dhdz != 0.0:
            DZ = -H4 / (cfg.ENMX2 * val_dhdz)
        else:
            DZ = 0.0 # Should not happen

        # Update Enthalpy
        H_new = H + val_dhdz * DZ

        # Overshoot handling (End of LV Region)
        if H_new >= cfg.HV:
            # Calculate partial step to exactly hit HV
            if val_dhdz != 0.0:
                DZ = (cfg.HV - H) / val_dhdz
            
            H = cfg.HV
            state.H[jj] = H
            state.Z[jj] = state.Z[jj - 1] + DZ
            state.WFV[jj] = 1.0
            
            if cfg.print_rows:
                _print_lv_row(jj, j_ref, state.Z[jj], cfg.TVAP, H, 1.0)
            
            # Update state index and exit
            state.ii = jj
            return

        # Normal Step
        H = H_new
        state.H[jj] = H
        state.Z[jj] = state.Z[jj - 1] + DZ
        
        wfv = (H - cfg.HL) / denom if denom != 0.0 else 0.0
        state.WFV[jj] = wfv

        if cfg.print_rows:
            _print_lv_row(jj, j_ref, state.Z[jj], cfg.TVAP, H, wfv)

        jj += 1
        j_ref += 1

    # If we hit max_steps
    state.ii = jj