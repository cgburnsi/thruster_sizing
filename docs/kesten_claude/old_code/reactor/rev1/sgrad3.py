import numpy as np
from unbar import UNBAR
from trapp1 import TRAPP

def KCF(A, B, C, D, E):
    return (0.61*A)/B * (C/(B*D))**-0.667 * (A/(E*C))**-0.41

# Equation I-3: Diffusion Coefficient of reactant gas in the catalyst particle at operating conditions
def DP3F(X, Y, Z):
    return 14.7*Y/Z * (X/492.)**1.823 * (1.0-np.exp(-0.0672*Z*492.0/(14.7*X)))

def SGRAD(TEMP, PRES, GG, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3):
    
    # Initialization
    WAF1                = 0.8                               # [-] Weighting Factor #1
    WAF2                = 1.0 - WAF1                        # [-] Weighting Factor #2
    
    # Initialization
    LTFLG               = 0                                 # [-] Flag for something around line 870
    LP1                 = 1                                 # [-] Loop Variable
    TPSP                = 0.0                               # [degR] Holds a previous value of the TPS
    T                   = TEMP                              # [degR] Interstitial Fluid Temperature
    P                   = PRES                              # [psia] Interstitial Fluid Pressure
    G                   = GG                                # [lb/ft^2-sec] Mass Flow Rate
    NPART               = 50                                # [-] Number of Steps in the Trapazoidal solver
    ALPH2               = ALPHA2                            # [(lb/ft^3)^1.6/sec] Preexponential Factor (See main.py function)
    ALPH3               = ALPHA3                            # [1/sec] Preexponetial factor in the rate (See main.py function)
    CI1, CI2, CI3, CI4  = C1, C2, C3, C4                    # [lb/ft^3] Interstitial Species Concentrations
    D03, D04            = DIF3, DIF4                        # [ft2/s] Diffusion Coefficients (Gas Phase at STP)
    RHO                 = CI1 + CI2 + CI3 + CI4             # [lb/ft^3] Interstitial Fluid Density
    MU                  = UNBAR('VISVST', T)                # [lb/ft-s] Interstitial Fluid Viscosity (maybe N2H4 only?? Not Sure)
    CF1                 = UNBAR('CFTBL1', T)                # [BTU/lb-degR] Specific heat of hydrogen
    CF2                 = UNBAR('CFTBL2', T)                # [BTU/lb-degR] Specific heat of nitrogen
    CF3                 = UNBAR('CFTBL3', T)                # [BTU/lb-degR] Specific heat of ammonia
    CF4                 = UNBAR('CFTBL4', T)                # [BTU/lb-degR] Specific heat of hydrazine    
    
    # Heater Transfer Coefficient (from Ref. 1)
    CFBAR = (CI1*CF1 + CI2*CF2 + CI3*CF3 + CI4*CF4) / RHO   # [BTU/lb-degR] Average Specific Heat of Interstitial Fluid
    HC    = 0.74* G * CFBAR * (G / (AP*MU))**-0.41          # [BTU/ft2-sec-degR] Heater Transfer Coefficient

    # Diffusion Coefficients at Operating Conditions
    DI3 = D03 * (14.7/P) * (T/492.0)**1.823
    DI4 = D04 * (14.7/P) * (T/492.0)**1.823 
    
    # Mass Transfer Coefficients
    KC3 = KCF(G, RHO, MU, DI3, AP)
    KC4 = KCF(G, RHO, MU, DI4, AP)
    
    # Step 1: Initial Guess for the reactant concentration at the surface of the particle (Cpi/2)    
    DP3 = DP3F(T, D03, P)
    CPS = CI3 / 2.0
    
    # Step2 : Calculate a slope at the surface of the particle
    CMCPN = CI3 - CPS
    DCPDX = (KC3 / DP3) * (CI3 - CPS)     # [?] Equation I-2: Slope of concentration profile at the particle surface (x=a)

    # Step 3: Calculate the temperature at the particle Surface
    H3 = UNBAR('H3TBL', T)                                  # [BTU/lb] Heat of Reaction for ammonia (Varies with each heat iteration loop)
    H4 = UNBAR('H4TBL', T)                                  # [BTU/lb] Heat of Reaction for hydrazine (Constant for each entry to routine)
    
    if LP1 == 1:
        TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC  # [degC] Equation I-4: Temperature at Particle Surface
    else:
        TPSPP = TPSP
        TPSP  = TPS
        TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC  # [degC] Equation I-4: Temperature at Particle Surface    
    if TPS < 0: TPS = 1.0
    
    # Step 4: Linear Extrapolation  to determine X0
    H3    = UNBAR('H3TBL', TPS)                                # [BTU/lb] Heat of Reaction for ammonia at current particle surface temperature
    H3P   = H3                                                # [BTU/lb] Save current value for the next loop iteration
    DP3   = DP3F(TPS, D03, P)                                 # [ft2/s] Diffusion Coefficient at TPS
    DP3P  = DP3                                              # [ft2/s] Save current value for the next loop iteration
    TMTPN = T-TPS                                           # [degR] ?? Some type of temperature value for the iteration loop
    GAMMA = BGM/TPS                                         # [-???] Not sure what this is
    BETA  = -CPS*H3*DP3/(KP*TPS)
    K0    = ALPH2*np.exp(-GAMMA)*CI1**EN3                   # [?] Reaction Rate Constant	
    X0    = A - CPS/DCPDX                                      # [ft] Equation I-5: Location where reactant concentration profile changes
    X0P   = X0                                              # [ft] Stores the value of X0 for next iteration loop

    if X0 < 0:
        X0 = 0.0
        XOA = 0.0
        CPS = CI3/(DP3/(A*KC3)+1.0)
        DCPDX = CI3/A
        TPS = T-(H4*KC4*CI4+H3*DP3*DCPDX)/HC	
        if TPS < 0: TPS = 1.0
        print('We have calculated a negative X0 during iteration.')
        RIESUM = TRAPP(XOA, 1.0, NPART, XOA, CPS, GAMMA, BETA, K0)
    else:
        XOA = X0/A
        RIESUM = TRAPP(XOA, 1.0, NPART, XOA, CPS, GAMMA, BETA, K0)
    
    # Calculate the new CPS value
    CPSP  = CPS
    CMCPO = CMCPN
    CPS   = CI3 - A*RIESUM/KC3

    if LTFLG == 0:                                  # GO TO 84        
        LTFLG = 0
        if CPS < 0:
            ...
        else:
            CPS = 0.0
            CPS = 0.2*CPS + 0.8*CPSP
            DCPDX = KC3/DP3P*(CI3-CPS)
            CMCPN = CI3-CPS
            H3 = H3P
            LP1 = LP1+1
            if LP1-25 < 0:
                WAF1 = WAF1 + 0.05
            else:
                ...
    else:                                       # GO TO 80
        if CPS-(0.25*CI3) > 0:
            CMCPN = CI3 - CPS
            # Calculate a new TP
            DCPDX = KC3/DP3*(CI3-CPS)
            GRAD = DCPDX*DP3
            TGRAD = HC*(T-TPS)
            TPSPP = TPSP
            TPSP = TPS
            TMTPO = TMTPN
            TPS = T - (H4*KC4*CI4 + H3*DP3*DCPDX) / HC
            if TPS < 0.0: TPS = 1.0
            H3 = UNBAR('H3TBL', TPS)
            DP3 = DP3F(TPS,D03,P)
            TMTPN = T - TPS
            GAMMA = BGM/TPS
            BETA = -CPS*H3*DP3/(KP*TPS)
            K0 = ALPH3*np.exp(-GAMMA) * CI1**EN3
            
            # TEST FOR 5% CONVERGENCE
            VALUE1 = np.abs(TMTPO-TMTPN) / TMTPN
            VALUE2 = np.abs(CMCPO-CMCPN) / CMCPN
            
            print(VALUE1, VALUE2)
        else:
            LTFLG = 1
    
    

    