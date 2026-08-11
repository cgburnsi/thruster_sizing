import numpy as np
import sys                         # <-- ADD THIS LINE
import logging
from logger_setup import LoggerSetup # <-- ADD THIS LINE

# Import all functions from our modules
from print_data import print_inputs, print_liquid_data, print_lqvp_data, print_vapor_data, print_vapor_debug
from slope1 import SLOPE
from unbar import UNBAR
from conc01 import CONC
from param01 import PARAM
from sgrad8 import SGRAD
from lqvp1 import LQVP

# In main5.py...
# (SimulationParameters class is unchanged)

def VAPOR(H_in, Z_in, II_in, TEMP_in, PRES_in,
          C1_in, C2_in, C3_in, C4_in,
          Z, DERIV, DHDZ, # Note: DERIV and DHDZ arrays are no longer needed
          FC, HF, HL, HV, ENMX3, G0, TVAP, Z0, DZ_in,
          AGM, BGM, ALPHA1, ALPHA2, ALPHA3, CGM,
          DIF3, DIF4, PRES, KP, 
          WM1, WM2, WM3, WM4, R, ZEND,
          EN1, EN2, EN3, GATZ0, max_steps=1000): # <-- GATZ0 added
    """
    Solves the vapor region of the reactor.
    [cite_start]This is a Python implementation of SUBROUTINE VAPOR [cite: 3487-3567].
    """
    print(' ------------  Vapor Region ------------')

    # --- 1. Initialize State from LQVP ---
    II = II_in
    Z[II] = Z_in
    TEMP = TEMP_in
    PRES = PRES_in
    H = H_in
    C1 = C1_in
    C2 = C2_in
    C3 = C3_in
    C4 = C4_in
    DZ = DZ_in # Use the last step size from LQVP as a starting guess
    
    # --- 2. Calculate initial mol fractions and FRAC3D ---
    SUM = C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4
    FRAC1 = (C1/WM1) / SUM
    FRAC2 = (C2/WM2) / SUM
    FRAC3 = (C3/WM3) / SUM
    FRAC4 = (C4/WM4) / SUM
    if (3. - FRAC1/FRAC2) == 0.0:
        FRAC3D = 0.0
    else:
        FRAC3D = (FRAC1/FRAC2 - 1.) / (3. - FRAC1/FRAC2)
    
    # Print the initial vapor state
    print_vapor_data(Z[II], TEMP, PRES, H, C1, C2, C3, C4, FRAC3D)
    
    # --- 3. Start Vapor Integration Loop ---
    while II < (max_steps - 1): # Safety margin for Z[II+1]
        
        # --- 4. Get current properties at Z[II] and TEMP ---
        AP = UNBAR('ZTBLAP', Z[II])
        A = UNBAR('ZTBLA', Z[II])
        DELA = UNBAR('ZTBLD', Z[II])
        
        # --- 5. Get heat of reactions ---
        H1 = UNBAR('SHTBL1', TEMP) # H2 (is 0.0)
        H2 = UNBAR('SHTBL2', TEMP) # N2 (is 0.0)
        H3 = UNBAR('H3TBL', TEMP) # NH3
        H4 = UNBAR('TBLH4', TEMP)  # N2H4

        # --- 6. Get NH3 catalytic rate (T3) ---
        # Always call SGRAD. It has its own internal checks for low temp.
        G, GMMA, K, BETA, DPA_NH3 = PARAM(
            TEMP, Z[II], 1, C3, H3, 1, # LVOP=1 for vapor/NH3
            Z0, G0, FC, GATZ0,
            AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1
        )
        GRAD, TGRAD = SGRAD(TEMP, PRES, G, C1, C2, C3, C4, 
                            DIF3, DIF4, A, AP, 
                            GMMA, KP, ALPHA2, ALPHA3, EN3)
        
        T3 = GRAD * AP # Catalytic NH3 rate [lb/ft^3-s]

        # --- 7. Get N2H4 catalytic rate (T4) ---
        # Get N2H4 gas-phase parameters
        G, _, _, _, DPA_N2H4 = PARAM(
            TEMP, Z[II], 1, C4, H4, 1, # LVOP=1
            Z0, G0, FC, GATZ0,
            AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1
        )
        VIS = UNBAR('VISVST', TEMP)
        RHO = (C1+C2+C3+C4)
        if RHO == 0: RHO = 1e-10 # prevent divide by zero
        
        # (Eq. 10)
        KC4 = (0.616 / RHO) * (VIS / (RHO * DPA_N2H4))**(-0.667) * (G / (AP * VIS))**(-0.41)
        T4 = KC4 * C4 * AP # Catalytic N2H4 rate [lb/ft^3-s]
        
        # --- 8. Get N2H4 homogeneous thermal rate (T2) ---
        RHOM = ALPHA3 * C4 * np.exp(-CGM / TEMP)
        T2 = RHOM * DELA # [lb/ft^3-s]

        # --- 9. Get species specific heats (Cp) ---
        CP1 = UNBAR('CFTBL1', TEMP)
        CP2 = UNBAR('CFTBL2', TEMP)
        CP3 = UNBAR('CFTBL3', TEMP)
        CP4 = UNBAR('CFTBL4', TEMP)
        CFBAR = (C1*CP1 + C2*CP2 + C3*CP3 + C4*CP4) / max(RHO, 1e-10)
        
        # --- 10. Calculate Gradients ---
        
        if Z[II] < Z0:
            F = FC
        else:
            F = 0.0
        
        # (Eq. 9) Heat transfer coefficient
        HC = 0.74 * (G/ (AP*VIS))**(-0.41) * (CFBAR * G) 
        
        # --- Enthalpy Gradient Terms (Eq. 4) ---
        H_N2H4_Term = -H4 * (T2 + T4)  # Heat from N2H4 decomp
        H_NH3_Term  = -H3 * T3         # Heat from NH3 decomp
        HC_Loss_Term = -HC * AP * (TEMP - TGRAD)
        F_Loss_Term  = -F * (H - HF)
        
        # --- AFTER ---
        # H_NH3_Term is removed because it's already accounted for *inside* the
        # HC_Loss_Term, which is calculated based on TGRAD from SGRAD.
        DHDZ_val = (1.0 / G) * (H_N2H4_Term + HC_Loss_Term + F_Loss_Term)        
        # --- FIX: Overwrite T4 with the complex SLOPE result ---
        # We just called SLOPE (DERIF). Now we must *use* that value.
        # This is the real catalytic N2H4 rate [lb/ft^3-s].
        T4 = DERIF * DPA_N2H4 * AP
        
        # Pressure Gradient (Eq. 14)
        DPDZ = -((1.0 - DELA) / DELA**3) * (1.75 + 150.0 * (1.0 - DELA) / (2.0 * A * G / VIS)) * (G**2 / (32.2 * A * RHO * 144.0))

        # Species Weight Fraction Gradients (Eq. 5, 6, 7)
        DW4DZ = (1.0/G) * (F - T2 - T4 - C4 * F / RHO)
        DW3DZ = (1.0/G) * (T2*WM3/WM4 + T4*WM3/WM4 - T3 - C3 * F / RHO)
        DW2DZ = (1.0/G) * (0.5*T2*WM2/WM4 + 0.5*T4*WM2/WM4 + 0.5*T3*WM2/WM3 - C2 * F / RHO)
        DW1DZ = (1.0/G) * (0.5*T2*WM1/WM4 + 0.5*T4*WM1/WM4 + 1.5*T3*WM1/WM3 - C1 * F / RHO)
        
        # --- 11. Calculate new step size (DZ) ---
        '''
        if DHDZ_val == 0.0:
            DZ = 0.01 # Prevent divide by zero if rate is zero
        else:
            # --- FIX: Handle negative DHDZ (temp peak) ---
            if DHDZ_val < 0:
                # We have passed the temperature peak.
                # Use a small, fixed positive step to continue.
                DZ = 1e-4 # Use a small fixed step of 0.0001 ft
            else:
                # Original logic for when temp is rising
                DZ = -H4 / (ENMX3 * DHDZ_val)
        '''
        DZ = 1.0e-5
        
        # --- 12. Check for overshoot ---
        if Z[II] + DZ > ZEND:
            DZ = ZEND - Z[II]
            
        # --- 13. Safety Check for Negative/Zero DZ ---
        # (This should not be triggered by DHDZ anymore, but good to keep)
        if DZ <= 0.0:
            print(f"--- ERROR: Negative or Zero DZ detected ({DZ:.2e}) ---")
            print_vapor_debug(Z[II], TEMP, DHDZ_val, DZ, T2, T3, T4, H3, H4, TGRAD, HC_Loss_Term, F_Loss_Term)
            print("--- ABORTING SIMULATION ---")
            break

        # --- 14. DEBUG PRINT ---
        print_vapor_debug(Z[II], TEMP, DHDZ_val, DZ, T2, T3, T4, H3, H4, TGRAD, HC_Loss_Term, F_Loss_Term)
        
        # --- 15. Take Euler step ---
        Z[II+1] = Z[II] + DZ
        H = H + DHDZ_val * DZ
        PRES = PRES + DPDZ * DZ
        
        # --- 15a. Calculate new Weight Fractions ---
        # (RHO is the *old* density from the start of this step)
        W1_new = (C1 / RHO) + DW1DZ * DZ
        W2_new = (C2 / RHO) + DW2DZ * DZ
        W3_new = (C3 / RHO) + DW3DZ * DZ
        W4_new = (C4 / RHO) + DW4DZ * DZ
        
        # --- 15b. Calculate new Temperature ---
        # (Use CFBAR (CPBAR_old) from the start of this step for a stable Euler step)
        if CFBAR == 0.0:
            print(f"--- ERROR: CFBAR is zero. Aborting. ---")
            break
        TEMP = TEMP + (DHDZ_val * DZ) / CFBAR
        
        # --- 15c. Update Concentrations for the *next* loop iteration ---
        
        # Ensure weight fractions are not negative
        W1_new = max(W1_new, 0.0)
        W2_new = max(W2_new, 0.0)
        W3_new = max(W3_new, 0.0)
        W4_new = max(W4_new, 0.0)
        
        # Normalize weight fractions to ensure they sum to 1.0
        # This prevents drift in the Euler integration
        W_sum = W1_new + W2_new + W3_new + W4_new
        if W_sum > 0.0:
            W1_new /= W_sum
            W2_new /= W_sum
            W3_new /= W_sum
            W4_new /= W_sum
            
        # Calculate new average molecular weight
        SUM_INV_WM = (W1_new/WM1 + W2_new/WM2 + W3_new/WM3 + W4_new/WM4)
        if SUM_INV_WM == 0.0:
            print(f"--- ERROR: Sum of inverse molecular weights is zero. Aborting. ---")
            # This can happen if all W_new are 0
            WM_avg_new = 1.0 # arbitrary to prevent crash
        else:
            WM_avg_new = 1.0 / SUM_INV_WM

        # Calculate new Density (RHO_new) using Ideal Gas Law
        # R = 10.73 [psia-ft3/lb-mol-degR]
        RHO_new = PRES * WM_avg_new / (R * TEMP)
        
        # Set the loop variables for the *next* iteration
        C1 = W1_new * RHO_new
        C2 = W2_new * RHO_new
        C3 = W3_new * RHO_new
        C4 = W4_new * RHO_new
        
        # RHO will be recalculated at the top of the next loop, but
        # we need to set it here for the mol fraction calculation below
        RHO = RHO_new
            
        #TEMP = TEMP + (DHDZ_val * DZ) / CPBAR_new
        
        # --- 16. Update mol fractions and FRAC3D ---
        SUM = C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4
        FRAC1 = (C1/WM1) / SUM
        FRAC2 = (C2/WM2) / SUM
        FRAC3 = (C3/WM3) / SUM
        FRAC4 = (C4/WM4) / SUM
        
        if (3. - FRAC1/FRAC2) == 0.0:
            FRAC3D = 0.0
        else:
            FRAC3D = (FRAC1/FRAC2 - 1.) / (3. - FRAC1/FRAC2)

        # --- 17. Print and Increment ---
        II += 1
        print_vapor_data(Z[II], TEMP, PRES, H, C1, C2, C3, C4, FRAC3D)
        
        # --- 18. Check for end ---
        if Z[II] >= ZEND:
            print("--- Vapor Region Solution Complete (End of Bed) ---")
            break
            
    return H, II, TEMP


if __name__ == '__main__':
    
    # --- 1. Set up the logger ---
    # This will duplicate all 'print' output to 'simulation.log'
    sys.stdout = LoggerSetup('simulation.log')

    try:
        # --- 2. Initialize all variables ---
        max_steps   = 1000
        Z           = np.zeros(max_steps)
        DERIV       = np.zeros(max_steps)
        DHDZ        = np.zeros(max_steps)
        DZ          = 0.0
        C1          = 1.0 # Note: Changed to float
        
        # Fluid Properties
        HF = 0.0
        TF = 530.0
        CFL = 0.7332
        WM1 = 2.016
        WM2 = 28.016
        WM3 = 17.032
        WM4 = 32.048
        R = 10.73
        DIF3 = 0.17e-3
        DIF4 = 0.95e-4
        ALPHA1 = 1.0e10
        ALPHA2 = 1.0e11
        ALPHA3 = 2.14e10
        EN1 = 1.0
        EN2 = 1.0
        EN3 = -1.6
        ENMX1 = 200.0 # Note: Changed to float
        ENMX2 = 40.0
        ENMX3 = 80.0

        # Reactor Variables
        G0 = 3.00
        Z0 = 0.0
        ZEND = 0.25
        FC = 0.0
        AGM = 2500.0
        BGM = 50000.0
        CGM = 33000.0
        PRES = 100.0
        KP = 0.4e-4

        # Initialization
        IFC = 1
        TVAP = UNBAR('TVAP', PRES)
        
        HL = 212.628
        HV = 715.478
        
        GATZ0 = G0 + FC * Z0
        H = HF
        
        if FC <= 0: IFC=0
        
        # =========================================================================
        # --- 3. SOLVE LIQUID REGION ---
        # =========================================================================
        print("\n--- Solving Liquid Region ---")
        
        II = 0
        Z[II] = 0.0
        DZ = 0.0
        
        while H < HL:
            if II > 0:
                Z[II] = Z[II-1] + DZ
                DHDZ[II] = DHDZ[II-1] 
                
            if Z[II] > ZEND:
                print(f"Error: Did not reach boiling point within bed length.")
                break
                
            TEMP = TF + (H - HF) / CFL
            VP = UNBAR('TBLVP', TEMP)
            CN2H4 = (VP * WM4) / (R * TEMP)
            H4 = UNBAR('TBLH4', TEMP)
            AP = UNBAR('ZTBLAP', Z[II])
            A = UNBAR('ZTBLA', Z[II])

            G, GMMA, K, BETA, DPA = PARAM(
                TEMP, Z[II], 1, CN2H4, H4, 0, Z0, G0, FC, GATZ0,
                AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1
            )
            
            DERIV[II], MI = SLOPE(CN2H4, GMMA, K, BETA, DPA, A, EN1, HL, DIF3, TEMP, PRES, G, AP, WM4, R)

            DHDZ[II] = -(H4 * DPA * AP * DERIV[II] + FC * (H - HF)) / G
            DZ = -H4 / (ENMX1 * DHDZ[II])
            
            print_liquid_data(Z[II], TEMP, H, DHDZ[II])

            if DHDZ[II] <= 0:
                print(f"Error: Integration stalled at Z={Z[II]:.4e}. DHDZ is not positive.")
                break

            if H + DHDZ[II] * DZ > HL:
                DZ = (HL - H) / DHDZ[II] 
                H = HL
            else:
                H = H + DHDZ[II] * DZ
                
            II += 1
            if II >= max_steps:
                print("Error: Reached max_steps before boiling.")
                break

        ZLV = Z[II-2] + DZ 
        Z[II-1] = ZLV
        TEMP = TF + (H - HF) / CFL
        DHDZ[II-1] = DHDZ[II-2] 
        print_liquid_data(Z[II-1], TEMP, H, DHDZ[II-1])
        
        print(f"\nLiquid Region Complete.")
        print(f"Entering LQVP at Z = {ZLV:.8e}, expecting ~7.37e-04 ft")
        
        # =========================================================================
        # --- 4. SOLVE LIQUID-VAPOR REGION ---
        # =========================================================================
        print(' ------------  Liquid-Vapor Region ------------')

        H_lv_end, II_lv_end, TEMP_lv_end = LQVP(
             H, ZLV, DERIV[II-1], II, DHDZ[II-1], TEMP, 
             DERIV, DHDZ, Z, 
             FC, HF, HL, HV, ENMX2, G0, TVAP, Z0, DZ, 
             AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1, 
             max_steps=250
        )
        
        # =========================================================================
        # --- 5. SOLVE VAPOR REGION ---
        # =========================================================================
        
        C1_init, C2_init, C3_init, C4_init = CONC(
            TEMP_lv_end, PRES, WM4, WM3, WM2, WM1, R, H_lv_end, HF
        )
        
        VAPOR(
            H_lv_end, Z[II_lv_end-1], II_lv_end, TEMP_lv_end, PRES,
            C1_init, C2_init, C3_init, C4_init,
            Z, DERIV, DHDZ, 
            FC, HF, HL, HV, ENMX3, G0, TVAP, Z0, DZ, 
            AGM, BGM, ALPHA1, ALPHA2, ALPHA3, CGM,
            DIF3, DIF4, PRES, KP, 
            WM1, WM2, WM3, WM4, R, ZEND,
            EN1, EN2, EN3, GATZ0, max_steps=max_steps
        )

    finally:
        # --- 6. SHUTDOWN LOGGER ---
        # This ensures the log file is saved, even if the script crashes
        logging.shutdown()
        # Restore the original stdout
        sys.stdout = sys.__stdout__




