import numpy as np

from print_data import print_inputs, print_liquid_data
from slope1 import SLOPE
from unbar import UNBAR
from conc01 import CONC
from param01 import PARAM
from sgrad8 import SGRAD
from lqvp1 import LQVP


if __name__ == '__main__':
    
    # Simulation Variables
    max_steps   = 1000                      # [-] Maximum Number of Steps Allowed in Calculation
    Z           = np.zeros(max_steps)       # [ft] Axial stations
    DERIV       = np.zeros(max_steps)
    DHDZ        = np.zeros(max_steps)       # [ft] Step length along reactor length.  Calculated each loop
    DZ          = 0.0                       # [ft] Length step to the next calculation point.  Calculated each loop
    C1          = 1                         # [?] Placeholder for vapor-phase logic
    
    # Fluid Properties
    HF = 0                  # [Btu/lb] Enthalpy of liquid hydrazine entering the bed
    TF = 530                # [degR] Temperature of liquid hydrazine entering the bed
    CFL = 0.7332            # [Btu/lb-degR] Specific Heat of Liquid Hydrazine
    WM1 = 2.016             # [lb/lb-mol] Molecular Weight of Hydrogen
    WM2 = 28.016            # [lb/lb-mol] Molecular Weight of Nitrogen
    WM3 = 17.032            # [lb/lb-mol] Molecular Weight of Ammonia
    WM4 = 32.048            # [lb/lb-mol] Molecular Weight of Hydrazine
    R = 10.73               # [psia-ft3/lb-mol-degR] Gas Constant
    DIF3 = 0.17e-3          # [ft2/s] Diffusion coefficient of ammonia in the gas phase at STP
    DIF4 = 0.95e-4          # [ft2/s] Diffusion coefficient of hydrazine in the gas phase at STP
    ALPHA1 = 1.0e10         # [1/sec] Preexponetial factor in the rate equation for the catalytic decomposition of hydrazine
    ALPHA2 = 1.0e11         # [(lb/ft^3)^1.6/sec] Preexponential factor in the rate equation for the catalytic decomposition of ammonia
    ALPHA3 = 2.14e10        # [1/sec] Preexponetial factor in the rate equation for the thermal decomposition of hydrazine
    EN1 = 1.0               # [-] Order of hydrazine catalytic decomposition reaction with to hydrazine hydrazine
    EN2 = 1.0               # [-] Order of ammonia catalytic decomposition reaction with respect to ammonia
    EN3 = -1.6              # [-] Order of ammonia catalytic decomposition reaction with respect to hydrogen
    ENMX1 = 200             # [-] Constant used to determine axial station increments in liquid region
    ENMX2 = 40              # [-] constant used to determine axial station increments in the liquid-vapor region
    ENMX3 = 80              # [-] constant used to determine axial station increments in the vapor region


    # Reactor Variables
    G0 = 3.00               # [lb/ft2-s] Inlet Mass Flow Rate
    Z0 = 0.0                # [ft] Axial Distance to the End of a Buried Injector
    ZEND = 0.25             # [ft] Bed Length
    FC = 0                  # [lb/ft3-sec] Rate of Feed of Hydrazine into System
    AGM = 2500              # [degR] Activation energy for catalytic decomposition of hydrazine divided by the gas constant
    BGM = 50000             # [degR] Activation energy for catalytic decomposition of ammonia divided by the gas constant
    CGM = 33000             # [degR] Activation energy for thermal decomposition of hydrazine divided by the gas constant
    PRES = 100.0            # [psia] Inlet Chamber Pressure
    KP = 0.4e-4             # [Btu/ft-sec-degR] Thermal Conductivity of the porous catalyst particle (Shell 405)


    
    
    # Initialization
    IFC = 1
    TVAP = UNBAR('TVAP', PRES)
    DELHV = UNBAR('DHVST', TVAP)
    DELHL = UNBAR('DHLVST', TVAP)
    #HL = (TVAP - TF) * CFL      # Enthalpy at saturated liquid
    #HV = HL + DELHV - DELHL     # Enthalpy at saturated vapor
    
    HL = 212.628              # [BTU/lb] From report output log (Fig 3c)
    HV = 715.478              # [BTU/lb] From report output log (Fig 3d)
    
    GATZ0 = G0 + FC * Z0
    H = HF                      # Start at feed enthalpy
    
    if FC <= 0: IFC=0
    
    #print_inputs(HF, HL, HV, TF, TVAP, CFL, PRES, KP, FC, G0, R, ALPHA3, CGM, DIF3, DIF4, WM4, WM3, WM2, WM1, ZEND, AGM, BGM, ALPHA1, ALPHA2, EN1, EN2, EN3, ENMX1, ENMX2, ENMX3, Z0)
    
    # Liquid Region
    # Main Loop (Replaces lines 67-130)
    print("\n--- Solving Liquid Region ---")
    
    II = 0 # Start at index 0
    Z[II] = 0.0
    DZ = 0.0 # DZ is calculated *after* the state is evaluated
    
    while H < HL:
        # --- Store current state (for all but the very first point) ---
        if II > 0:
            Z[II] = Z[II-1] + DZ
            DHDZ[II] = DHDZ[II-1] # Store the DHDZ that *led* to this state
            
        # --- Check for Bed End ---
        if Z[II] > ZEND:
            print(f"Error: Did not reach boiling point within bed length.")
            print(f"Reached Z = {Z[II]:.4e} ft")
            break
            
        # 1. Get current properties at Z[II] and TEMP
        TEMP = TF + (H - HF) / CFL # (Eq. [cite_start]2) [cite: 98, 3420]
        VP = UNBAR('TBLVP', TEMP)
        CN2H4 = (VP * WM4) / (R * TEMP)
        H4 = UNBAR('TBLH4', TEMP)
        AP = UNBAR('ZTBLAP', Z[II])
        A = UNBAR('ZTBLA', Z[II])

        # 2. Get parameters
        G, GMMA, K, BETA, DPA = PARAM(
            TEMP, Z[II], 1, CN2H4, H4, 0, Z0, G0, FC, GATZ0,
            AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1
        )
        
        # [cite_start]3. Call SLOPE to get the reaction rate (DERIV) [cite: 3418-3423]
        DERIV[II], MI = SLOPE(CN2H4, GMMA, K, BETA, DPA, A, EN1, HL, DIF3, TEMP, PRES, G, AP, WM4, R)

        # 4. Calculate the enthalpy gradient (DHDZ) and the next step size (DZ)
        DHDZ[II] = -(H4 * DPA * AP * DERIV[II] + FC * (H - HF)) / G
        DZ = -H4 / (ENMX1 * DHDZ[II]) # [cite: 3413-3414]
        
        # [cite_start]5. Print output for the current step [cite: 3421]
        print_liquid_data(Z[II], TEMP, H, DHDZ[II])

        # 6. Check for stall (DHDZ is not positive)
        if DHDZ[II] <= 0:
            print(f"Error: Integration stalled at Z={Z[II]:.4e}. DHDZ is not positive.")
            break

        # 7. Check for overshoot (don't step past the boiling point)
        if H + DHDZ[II] * DZ > HL:
            DZ = (HL - H) / DHDZ[II] # Calculate partial step to land *exactly* on HL
            H = HL
        else:
            H = H + DHDZ[II] * DZ # Take the full step
            
        # 8. Increment index and check loop limit
        II += 1
        if II >= max_steps:
            print("Error: Reached max_steps before boiling.")
            break
    # --- End of Loop ---
    ZLV = Z[II-2] + DZ # Z[II-2] was the last point, DZ was the final step
    Z[II-1] = ZLV
    TEMP = TF + (H - HF) / CFL # This is now TVAP
    DHDZ[II-1] = DHDZ[II-2] # Carry over last DHDZ
    print_liquid_data(Z[II-1], TEMP, H, DHDZ[II-1]) # Print the final state
    
    print(f"\nLiquid Region Complete.")
    print(f"Entering LQVP at Z = {ZLV:.8e}, expecting ~7.37e-04 ft")
    
    # Liquid-Vapor Region
    print(' ------------  Liquid-Vapor Region ------------')

    # --- THIS IS THE FIX ---
    # Pass the LAST calculated values (at index II-1)
    # Pass the NEXT index (II) to start the LQVP loop
    LQVP(H, ZLV, DERIV[II-1], II, DHDZ[II-1], TEMP, 
         DERIV, DHDZ, Z, 
         FC, HF, HL, HV, ENMX2, G0, TVAP, Z0, DZ, 
         AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1, 
         max_steps=250)