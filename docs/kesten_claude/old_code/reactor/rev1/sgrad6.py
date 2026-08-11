import numpy as np
from unbar import UNBAR
from trapp1 import TRAPP

def KCF(A, B, C, D, E):
    return (0.61*A)/B * (C/(B*D))**-0.667 * (A/(E*C))**-0.41

# Equation I-3: Diffusion Coefficient of reactant gas in the catalyst particle at operating conditions
def DP3F(X, Y, Z):
    return 14.7*Y/Z * (X/492.)**1.823 * (1.0-np.exp(-0.0672*Z*492.0/(14.7*X)))

def step1(CI3):
    return CI3 / 2.0
    
def step2(T, D03, D04, P, CI3, CPS, G, RHO, MU, AP):
    
    # Diffusion Coefficients at Operating Conditions
    DI3 = D03 * (14.7/P) * (T/492.0)**1.823
    DI4 = D04 * (14.7/P) * (T/492.0)**1.823 
    
    # Mass Transfer Coefficients
    KC3 = KCF(G, RHO, MU, DI3, AP)
    KC4 = KCF(G, RHO, MU, DI4, AP)

    DP3 = DP3F(T, D03, P)
    DCPDX = (KC3 / DP3) * (CI3 - CPS)     # [?] Equation I-2: Slope of concentration profile at the particle surface (x=a)
    
    return DCPDX, DP3, KC3, KC4




def step3(T, G, MU, AP, CI1, CI2, CI3, CI4, RHO, CF1, CF2, CF3, CF4, KC4, DP3, DCPDX):
    
    H3 = UNBAR('H3TBL', T)                                  # [BTU/lb] Heat of Reaction for ammonia (Varies with each heat iteration loop)
    H4 = UNBAR('H4TBL', T)                                  # [BTU/lb] Heat of Reaction for hydrazine (Constant for each entry to routine)
    
    # Heater Transfer Coefficient (from Ref. 1)
    CFBAR = (CI1*CF1 + CI2*CF2 + CI3*CF3 + CI4*CF4) / RHO   # [BTU/lb-degR] Average Specific Heat of Interstitial Fluid
    HC    = 0.74* G * CFBAR * (G / (AP*MU))**-0.41          # [BTU/ft2-sec-degR] Heater Transfer Coefficient

    TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC  # [degC] Equation I-4: Temperature at Particle Surface
    return TPS
 



   
def SGRAD(TEMP, PRES, GG, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3):
    
    T                   = TEMP                              # [degR] Interstitial Fluid Temperature
    P                   = PRES                              # [psia] Interstitial Fluid Pressure
    G                   = GG                                # [lb/ft^2-sec] Mass Flow Rate
    D03, D04            = DIF3, DIF4                        # [ft2/s] Diffusion Coefficients (Gas Phase at STP)
    CI1, CI2, CI3, CI4  = C1, C2, C3, C4                    # [lb/ft^3] Interstitial Species Concentrations
    RHO                 = CI1 + CI2 + CI3 + CI4             # [lb/ft^3] Interstitial Fluid Density
    MU                  = UNBAR('VISVST', T)                # [lb/ft-s] Interstitial Fluid Viscosity (maybe N2H4 only?? Not Sure)
    CF1                 = UNBAR('CFTBL1', T)                # [BTU/lb-degR] Specific heat of hydrogen
    CF2                 = UNBAR('CFTBL2', T)                # [BTU/lb-degR] Specific heat of nitrogen
    CF3                 = UNBAR('CFTBL3', T)                # [BTU/lb-degR] Specific heat of ammonia
    CF4                 = UNBAR('CFTBL4', T)                # [BTU/lb-degR] Specific heat of hydrazine  


    # Step 1: Guess reactant concentration at the surface of the catalyst particle
    CPS = step1(CI3)
    
    # Step 2: Calculate a slope at the surface of the particle    
    DCPDX, DP3, KC3, KC4 = step2(T, D03, D04, P, CI3, CPS, G, RHO, MU, AP)
    
    CMCPN = CI3 - CPS

    # Step 3: Calculate the temperature at the particle Surface, TPS
    TPS = step3(T, G, MU, AP, CI1, CI2, CI3, CI4, RHO, CF1, CF2, CF3, CF4, KC4, DP3, DCPDX)
    
    # Step 4: Linear Extrapolation  to determine X0
    # Step I-5: Calculate X0
    # Step I-6: Numerical Integration (Trapezoidal Rule) of equation I-7 for new value of CPS
    # Step I-7: Calculate new DCPDX (equation I-3) using the new CPS caluclated in Step I-6
    # Step I-8: Calculate new TPS, DP3, GAMMA, BETA, and K0
    # Step I-9: Check for convergence

    
    
    
    
    
    
    # Phase II
    
    
    
    
    
    return 0
    
    
