import numpy as np
from unbar import UNBAR
from trapp1 import TRAPP

def KCF(A, B, C, D, E):
    return (0.61*A)/B * (C/(B*D))**-0.667 * (A/(E*C))**-0.41

# Equation I-3: Diffusion Coefficient of reactant gas in the catalyst particle at operating conditions
def DP3F(X, Y, Z):
    return 14.7*Y/Z * (X/492.)**1.823 * (1.0-np.exp(-0.0672*Z*492.0/(14.7*X)))

def CALC_KC(T, P, D03, D04, G, RHO, MU, AP):

    # Diffusion Coefficients at Operating Conditions
    DI3 = D03 * (14.7/P) * (T/492.0)**1.823
    DI4 = D04 * (14.7/P) * (T/492.0)**1.823 
    
    # Mass Transfer Coefficients
    KC3 = KCF(G, RHO, MU, DI3, AP)
    KC4 = KCF(G, RHO, MU, DI4, AP)    
    
    return KC3, KC4
    
def CALC_HC(T, G, AP, CI1, CI2, CI3, CI4):
    MU                  = UNBAR('VISVST', T)                # [lb/ft-s] Interstitial Fluid Viscosity (maybe N2H4 only?? Not Sure)
    CF1                 = UNBAR('CFTBL1', T)                # [BTU/lb-degR] Specific heat of hydrogen
    CF2                 = UNBAR('CFTBL2', T)                # [BTU/lb-degR] Specific heat of nitrogen
    CF3                 = UNBAR('CFTBL3', T)                # [BTU/lb-degR] Specific heat of ammonia
    CF4                 = UNBAR('CFTBL4', T)                # [BTU/lb-degR] Specific heat of hydrazine  
    RHO                 = CI1 + CI2 + CI3 + CI4             # [lb/ft^3] Interstitial Fluid Density

    # Heater Transfer Coefficient (from Ref. 1)
    CFBAR = (CI1*CF1 + CI2*CF2 + CI3*CF3 + CI4*CF4) / RHO   # [BTU/lb-degR] Average Specific Heat of Interstitial Fluid
    HC    = 0.74* G * CFBAR * (G / (AP*MU))**-0.41          # [BTU/ft2-sec-degR] Heater Transfer Coefficient
    
    return HC, RHO, MU



def CALC_TPS(T, KC3, KC4, CI4, DP3, HC, DCPDX):
    
    H3 = UNBAR('H3TBL', T)                                  # [BTU/lb] Heat of Reaction for ammonia (Varies with each heat iteration loop)
    H4 = UNBAR('H4TBL', T)                                  # [BTU/lb] Heat of Reaction for hydrazine (Constant for each entry to routine)

    TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC  # [degC] Equation I-4: Temperature at Particle Surface
    if TPS <= 0.0:
        print(f'Negative Surface Temperature Calculate TPS = {TPS:.2f} [degR]. TPS reset to 1.0')
        TPS = 1.0
    return TPS
 



   
def SGRAD(TEMP, PRES, GG, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3):    
    MAXITER             = 1                               # [-] Specify the maximum number of iterations to allow
    T                   = np.zeros(MAXITER)
    CPS                 = np.zeros(MAXITER)
    KC3                 = np.zeros(MAXITER)
    KC4                 = np.zeros(MAXITER)
    DP3                 = np.zeros(MAXITER)
    TPS                 = np.zeros(MAXITER)
    H3                  = np.zeros(MAXITER)
    X0                  = np.zeros(MAXITER)
    X0A                 = np.zeros(MAXITER)
    CMCPN               = np.zeros(MAXITER)
    CMCPO               = np.zeros(MAXITER)
    TMTPN               = np.zeros(MAXITER)
    TMTPO               = np.zeros(MAXITER)

    
    T[0]                = TEMP                              # [degR] Interstitial Fluid Temperature
    P                   = PRES                              # [psia] Interstitial Fluid Pressure
    G                   = GG                                # [lb/ft^2-sec] Mass Flow Rate
    NPART               = 50                                # [-] Number of Steps in the Trapazoidal solver
    ALPH2               = ALPHA2                            # [(lb/ft^3)^1.6/sec] Preexponential Factor (See main.py function)
    #ALPH3               = ALPHA3                            # [1/sec] Preexponetial factor in the rate (See main.py function)
    D03, D04            = DIF3, DIF4                        # [ft2/s] Diffusion Coefficients (Gas Phase at STP)
    CI1, CI2, CI3, CI4  = C1, C2, C3, C4                    # [lb/ft^3] Interstitial Species Concentrations
    LTFLG               = 0
    WAF1                = 0.8                               # Weight Factor
    WAF2                = 1 - WAF1
    TPSP                = 0
    
    CPS[0] = CI3 / 2.0                              # Step 1: Guess Reactant Concentration at the Surface of Particle, CPS
    
    for i in range(0, MAXITER):
        
        # Step 2: Calculate the concentratoin slope at the surface of the particle 
        HC, RHO, MU     = CALC_HC(T[i], G, AP, CI1, CI2, CI3, CI4)      # Calculate heat transfer coefficient, density, and viscosity 
        KC3[i], KC4[i]  = CALC_KC(T[i], P, D03, D04, G, RHO, MU, AP)    # Calculate the mass transfer coefficients at the operating conditions
        DP3[i]          = DP3F(T[i], D03, P)                            # Calculate ammonia diffusion coefficient
        DCPDX           = (KC3[i] / DP3[i]) * (CI3 - CPS[i])            # [?] Equation I-2: Slope of concentration profile at the particle surface (x=a)

        CMCPN[i] = CI3 - CPS[i]
        
        # Step 3: Calculate the temperature at the particle surface, TPS
        TPS[i]          = CALC_TPS(T[i], KC3[i], KC4[i], CI4, DP3[i], HC, DCPDX)    # [degR] Calculate the temperature at particle surface
    
        # Step 4: Linear Extrapolation  to determine X0
        H3[i]           = UNBAR('H3TBL', TPS[i])                   # [BTU/lb] Heat of Reaction for ammonia at current TPS
        DP3[i]          = DP3F(TPS[i], D03, P)                  # [ft2/s] Update diffusion coefficient
        H3P = H3[i]
        DP3P = DP3[i]
        
        TMTPN[i] = T[i]-TPS[i]                                           # [degR] ?? Some type of temperature value for the iteration loop
        GAMMA = BGM/TPS[i]                                         # [-???] Not sure what this is
        BETA  = -CPS[i]*H3[i]*DP3[i]/(KP*TPS[i])
        KO    = ALPH2*np.exp(-GAMMA)*CI1**EN3                   # [?] Reaction Rate Constant	

        # Step 5: Calculate X0
        X0P      = X0[i]
        X0[i]    = A - CPS[i]/DCPDX                                # [ft] Equation I-5: Location where reactant concentration profile changes
        X0A[i]   = X0[i] / A                                    # [?] oh boy.  this is not in the code.  they actually don't have the X0A initalized anywhere, so I'm going to take a guess and do this.
        if X0[i] < 0.0:
            print(f'Negative Location for Concentration Slope Change, X0 = {X0:.2f}, X0 reset to 0.0')
            X0[i]   = 0.0
            X0A[i]  = 0.0
            CPS[i]  = CI3/(DP3[i]/(A*KC3[i])+1.0)
            DCPDX   = CI3/A
            TPS[i]  = CALC_TPS(T[i], KC3[i], KC4[i], CI4, DP3[i], HC, DCPDX)    # [degR] Calculate the temperature at particle surface
		
        # Step 6: Calculate a new value for CPS (use trapezoidal rule)
        RIESUM = TRAPP(X0A[i], 1.0, NPART, X0A[i], CPS, GAMMA, BETA, KO)
        CPSP = CPS[i]
        CMCPO[i] = CMCPN[i]
        print(CPS)
        CPS[i] = CI3 - A*RIESUM / KC3[i]
        print(CPS)
        # Check CPS threshold and update LTFLG based on concentration
        THOLD = 0.25 * CI3
        if LTFLG == 0:
            if CPS[i] < THOLD:
                LTFLG = 1
                XOO = WAF1 * X0P + WAF2 * X0[i]
                CPS[i] = CI3 / (1.0+DP3[i]/(KC3[i]*A-KC3[i]*XOO))
                DCPDX = KC3[i] / DP3P * (CI3-CPS[i])
                CMCPN = CI3 - CPS[i]
                H3[i] = H3P
                continue
            else:
                print('we are good to go #1')               
        elif LTFLG == 1:
            LTFLG = 0
            if CPS[i] < 0.0:
                CPS[i] = 0.0
                CPS[i] = .2*CPS[i] + 0.8*CPSP
                DCPDX = KC3[i] / DP3P * (CI3-CPS[i])
                CMCPN = CI3 - CPS[i]
                H3[i] = H3P
                continue
            else:
                print('we are good to go. #2')      

        if i == 25:
            WAF1 += 0.05
            WAF2 = 1-WAF1
            if WAF1 >= 1.0:
                print(f'No Convergence with present weight factors WAF1 = {WAF1:.2f}, WAF2 = {WAF2:.2f}.  Adjusting weight factors')
        
    CMCPN[i] = CI3 - CPS[i]
    # Calculate a new TPS

    DCPDX = KC3[i]/DP3[i]*(CI3-CPS[i])
    GRAD = DCPDX * DP3[i]
    TGRAD = HC * (T-TPS)
    TPSPP = TPSP
    TPSP = TPS
    TMTPO[i] = TMTPN[i]
    TPS[i]  = CALC_TPS(T[i], KC3[i], KC4[i], CI4, DP3[i], HC, DCPDX)    # [degR] Calculate the temperature at particle surface
    H3[i] = UNBAR('H3TBL', TPS[i])
    DP3[i] = DP3F(TPS[i], D03, P)
    TMTPN[i] = T[i] - TPS[i]
    GAMMA = BGM / TPS[i]
    BETA = -CPS[i] * H3[i] * DP3[i] / (KP * TPS[i])
    KO    = ALPH2*np.exp(-GAMMA)*CI1**EN3                   # [?] Reaction Rate Constant	

    # Step 9: Check for Convergence
    CRIT1 = np.abs(TMTPO[i] - TMTPN[i])/TMTPN[i]
    #CRIT2 = np.abs(CMCPO[i] - CMCPN[i])/CMCPN[i]
    
    print(CMCPO[i], CMCPN[i], CI3, CPS[i])

    
    #print(CRIT1, CRIT2)
    print('we are good to go for the rest of things.')
                
            
            
        
        
    #print(RIESUM, GAMMA, BETA, CPS[i], X0[i], X0A[i])
            
    return 
    

    
    
    
    
    
    
