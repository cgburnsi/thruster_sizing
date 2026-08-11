import numpy as np

from print_data import print_inputs, print_liquid_data
from slope8 import SLOPE
from unbar import UNBAR
#from conc01 import CONC
from param01 import PARAM



if __name__ == '__main__':
    
    # Simulation Variables
    max_steps   = 2                       # [-] Maximum Number of Steps Allowed in Calculation
    Z           = np.zeros(max_steps)       # [ft] Axial stations
    DERIV       = np.zeros(max_steps)
    DHDZ        = np.zeros(max_steps)       # [ft] Step length along reactor length.  Calculated each loop
    DZ          = 0.0                       # [ft] Length step to the next calculation point.  Calculated each loop
    C1          = 1                         # [?] I'm not sure what this is.  It's part of PARAM, but I'm guessing it is determined later in the L-V phase calculation and not in the Liquid phase.
    
    # Fluid Properties
    HF = 0                  # [Btu/lb] Enthalpy of liquid hydrazine entering the bed
    TF = 530                # [degR] Temperature of liquid hydrazine entering the bed
    CFL = 0.7332            # [Btu/lb-degR] Specific Heat of Liquid Hydrazine
    WM1 = 2.016             # [lb/lb-mol] Molecular Weight of Hydrogen
    WM2 = 28.016            # [lb/lb-mol] Molecular Weight of Nitrogen
    WM3 = 17.032            # [lb/lb-mol] Molecular Weight of Amonia
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
    HL = (TVAP - TF) * CFL
    HV = HL + DELHV - DELHL
    GATZ0 = G0 + FC * Z0
    H = HF
    
    if FC <= 0: IFC=0
    
    #print_inputs(HF, HL, HV, TF, TVAP, CFL, PRES, KP, FC, G0, R, ALPHA3, CGM, DIF3, DIF4, WM4, WM3, WM2, WM1, ZEND, AGM, BGM, ALPHA1, ALPHA2, EN1, EN2, EN3, ENMX1, ENMX2, ENMX3, Z0)
    
    # Main Loop
    II = 1
    while II < max_steps:
        Z[II] = Z[II-1] + DZ       
        TEMP = TF + (H - HF) / CFL # Update the temperature of the fluid based on the ratio of enthalpy/specific heat
        VP = UNBAR('TBLVP', TEMP)
        CN2H4 = (VP * WM4) / (R * TEMP)
        H4 = UNBAR('TBLH4',TEMP)
        AP = UNBAR('ZTBLAP',Z[II])
        A = UNBAR('ZTBLA',Z[II])
                
        G, GMMA, K, BETA, DPA = PARAM(TEMP, Z[II], 1, CN2H4, H4, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1)

        DERIV[II], CPA, XOA, MI = SLOPE(CN2H4, GMMA, K, BETA, EN1, A)
        #if (H-HL) == 0.0:
        #    if(MI >= 20): DERIV[II] = DERIV[II-1]
        #DZ = 0.01  # temporary until full DZ logic is restored
    
        #DHDZ[II] = -(H4*DPA*AP*DERIV[II] + FC*(H-HL))/G
        #DZ = -H4/(ENMX1*DHDZ[II]) 
        
        
        print_liquid_data(Z[II], TEMP, H, DHDZ[II])


        # Increment Counter
        II += 1
        
        # Calculate the next DZ
    
    
    

    
  