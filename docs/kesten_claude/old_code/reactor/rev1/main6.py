from print_data import print_vapor_data_1

import numpy as np
from liquid1 import LIQUID
from liquidvapor1 import LIQUIDVAPOR
from param2 import PARAM
from conc01 import CONC
#from vapor1 import VAPOR
from unbar import UNBAR
from dataclasses import dataclass, field

@dataclass
class SimConfig:
    # Fluid Properties
    WM1         = 2.016             # [lb/lb-mol] Molecular Weight of Hydrogen
    WM2         = 28.016            # [lb/lb-mol] Molecular Weight of Nitrogen
    WM3         = 17.032            # [lb/lb-mol] Molecular Weight of Ammonia
    WM4         = 32.048            # [lb/lb-mol] Molecular Weight of Hydrazine
    ALPHA1      = 1.0e10            # [1/sec] Preexponetial factor in the rate equation for the catalytic decomposition of hydrazine
    ALPHA2      = 1.0e11            # [(lb/ft^3)^1.6/sec] Preexponential factor in the rate equation for the catalytic decomposition of ammonia
    ALPHA3      = 2.14e10           # [1/sec] Preexponetial factor in the rate equation for the thermal decomposition of hydrazine
    EN1         = 1.0               # [-] Order of hydrazine catalytic decomposition reaction with to hydrazine hydrazine
    EN2         = 1.0               # [-] Order of ammonia catalytic decomposition reaction with respect to ammonia
    EN3         = -1.6              # [-] Order of ammonia catalytic decomposition reaction with respect to hydrogen
    AGM         = 2500              # [degR] Activation energy for catalytic decomposition of hydrazine divided by the gas constant
    BGM         = 50000             # [degR] Activation energy for catalytic decomposition of ammonia divided by the gas constant
    CGM         = 33000             # [degR] Activation energy for thermal decomposition of hydrazine divided by the gas constant
    CFL         = 0.7332            # [Btu/lb-degR] Specific Heat of Liquid Hydrazine
    R           = 10.73             # [psia-ft3/lb-mol-degR] Gas Constant
    DIF3        = 0.17e-3           # [ft2/s] Diffusion coefficient of ammonia in the gas phase at STP
    DIF4        = 0.95e-4           # [ft2/s] Diffusion coefficient of hydrazine in the gas phase at STP
    KP          = 0.4e-4            # [Btu/ft-sec-degR] Thermal Conductivity of the porous catalyst particle (Shell 405)
    ENMX1       = 200               # [-] Constant used to determine axial station increments in liquid region
    ENMX2       = 40                # [-] constant used to determine axial station increments in the liquid-vapor region
    ENMX3       = 80                # [-] constant used to determine axial station increments in the vapor region
    C1          = 1                 # [?] I'm not sure what this is.  It's part of PARAM, but I'm guessing it is determined later in the L-V phase calculation and not in the Liquid phase.

    # Reactor Properties
    FC          = 0                 # [lb/ft3-sec] Rate of Feed of Hydrazine into System
    G0          = 3.00              # [lb/ft2-s] Inlet Mass Flow Rate
    Z0          = 0.0               # [ft] Axial Distance to the End of a Buried Injector
    ZEND        = 0.25              # [ft] Bed Length
    HF          = 0                 # [Btu/lb] Enthalpy of liquid hydrazine entering the bed
    TF          = 530               # [degR] Temperature of liquid hydrazine entering the bed
    PRES        = 100.0             # [psia] Inlet Chamber Pressure

    # Simulation Properties
    LIQ_MAX    = 100                # [-] Maximum Number of steps in the Liquid Zone
    LV_MAX     = 100                # [-] Maximum Number of steps in the Liquid-Vapor Zone
    VAP_MAX    = 100                # [-] Maximum Number of steps in the Vapor Zone
    
@dataclass
class SimState:
    max_iter:int        = field(init=False)
    Z:np.ndarray        = field(init=False)
    H:np.ndarray        = field(init=False)
    TEMP:np.ndarray     = field(init=False)
    PRES:np.ndarray     = field(init=False)
    DHDZ:np.ndarray     = field(init=False)
    
    def __post_init__(self):
        self.max_iter   = 300
        self.Z          = np.zeros(self.max_iter)
        self.H          = np.zeros(self.max_iter)
        self.TEMP       = np.zeros(self.max_iter)
        self.PRES       = np.zeros(self.max_iter)
        self.DHDZ       = np.zeros(self.max_iter)
        self.DERIV      = np.zeros(self.max_iter)
        self.Z1         = np.zeros(self.max_iter)
        self.H1         = np.zeros(self.max_iter)
        self.TEMP1      = np.zeros(self.max_iter)
        self.PRES1      = np.zeros(self.max_iter)
        self.DHDZ1      = np.zeros(self.max_iter)
        self.T1 = np.zeros(self.max_iter)
        self.T2 = np.zeros(self.max_iter)
        self.P2 = np.zeros(self.max_iter)
        self.H2 = np.zeros(self.max_iter)
        self.Z2 = np.zeros(self.max_iter)
        self.DHDZ2 = np.zeros(self.max_iter)
        self.DERIV2 = np.zeros(self.max_iter)
        self.T3 = np.zeros(self.max_iter)
        self.T3 = np.zeros(self.max_iter)
        self.P3 = np.zeros(self.max_iter)
        self.H3 = np.zeros(self.max_iter)
        self.Z3 = np.zeros(self.max_iter)
        self.DHDZ3 = np.zeros(self.max_iter)
        self.DERIV3 = np.zeros(self.max_iter)



# -------------------  The following goes to the SGRAD Module Eventually -------------------------
def RHETF(A, B, C, D, E, N):
    return E * A**(1 - N) * B**N * np.exp(C * D * (1.0 - B/A) / (1.0 + D * (1 - B/A)))

def FOXI1(X, R):
    return X**2 * R

def CPXF(X, Y, Z):
    return (X - Y) / (1 - Y) * Z

def TRAPP(U, V, NPART, XOA, CPS, GAMMA, BETA, K0):
    N = NPART - 1
    PART = NPART
    H = (V - U) / PART
    SUM = 0.0

    #print("\n🔍 Entering TRAPP integration loop...")
    #print(f"{'i':>3} {'X':>10} {'CPX':>10} {'RHET':>10} {'FOXI1':>10}")

    # First point
    CPX1 = CPXF(U, XOA, CPS)
    RHET1 = RHETF(CPS, CPX1, GAMMA, BETA, K0, 1)
    TRM1 = FOXI1(U, RHET1) / 2.0
    #print(f"{0:>3} {U:10.5e} {CPX1:10.5e} {RHET1:10.5e} {FOXI1(U, RHET1):10.5e}")

    # Internal points
    xs = np.linspace(U + H, V - H, N - 1)
    for i, x in enumerate(xs, start=1):
        CPX = CPXF(x, XOA, CPS)
        RHET = RHETF(CPS, CPX, GAMMA, BETA, K0, 1)
        FX = FOXI1(x, RHET)
        #print(f"{i:>3} {x:10.5e} {CPX:10.5e} {RHET:10.5e} {FX:10.5e}")
        SUM += FX

    # Last point
    CPX2 = CPXF(V, XOA, CPS)
    RHET2 = RHETF(CPS, CPX2, GAMMA, BETA, K0, 1)
    TRM2 = FOXI1(V, RHET2) / 2.0
    #print(f"{N:>3} {V:10.5e} {CPX2:10.5e} {RHET2:10.5e} {FOXI1(V, RHET2):10.5e}")

    RIESUM = H * (TRM1 + SUM + TRM2)

    #print(f"\n📏 TRAPP complete. H = {H:.5e}, RIESUM = {RIESUM:.5e}\n")
    return RIESUM

def KCF(A, B, C, D, E):
    return (0.61 * A) / B * (C / (B * D)) ** -0.667 * (A / (E * C)) ** -0.41

def DP3F(X, Y, Z):
    return 14.7 * Y / Z * (X / 492.) ** 1.823 * (1.0 - np.exp(-0.0672 * Z * 492.0 / (14.7 * X)))

def CALC_KC(T, P, D03, D04, G, RHO, MU, AP):
    DI3 = D03 * (14.7 / P) * (T / 492.0) ** 1.823
    DI4 = D04 * (14.7 / P) * (T / 492.0) ** 1.823
    KC3 = KCF(G, RHO, MU, DI3, AP)
    KC4 = KCF(G, RHO, MU, DI4, AP)
    return KC3, KC4

def CALC_HC(T, G, AP, CI1, CI2, CI3, CI4):
    MU = UNBAR('VISVST', T)
    CF1 = UNBAR('CFTBL1', T)
    CF2 = UNBAR('CFTBL2', T)
    CF3 = UNBAR('CFTBL3', T)
    CF4 = UNBAR('CFTBL4', T)
    RHO = CI1 + CI2 + CI3 + CI4
    CFBAR = (CI1 * CF1 + CI2 * CF2 + CI3 * CF3 + CI4 * CF4) / RHO
    HC = 0.74 * G * CFBAR * (G / (AP * MU)) ** -0.41
    return HC, RHO, MU


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

    if LTFLG != 1:
        if CPS <= 0.25 * CI3:
            LTFLG = 1
            X00 = WAF1*X0P+WAF2*X0
            CPS = CI3/(1.+DP3/(KC3*A-KC3*X00))
            DCPDX = KC3/DP3P*(CI3-CPS)
            CMCPN = CI3-CPS
            H3 = H3P
            LP1 = LP1+1
        else:
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
            print('There is more stuff here for lines 1100, 1110, 1120, and 1130')
    else:
        LTFLG = 0
        
        
        
        

    

    return GRAD, TGRAD






if __name__ == '__main__':
    
    # Initialization
    config  = SimConfig()
    state   = SimState()
    
    print('----- LIQUID ZONE -------------')
    state.Z1, state.T1, state.P1, state.H1, state.DHDZ1, state.DERIV1 = LIQUID(state, config)    
    print('----- LIQUID - VAPOR ZONE -----')
    state.Z2, state.T2, state.P2, state.H2, state.DHDZ2, state.DERIV2 = LIQUIDVAPOR(state, config)
    print('----- VAPOR ZONE --------------')
    
    
    # Retrieve Config Parameters
    Z0      = config.Z0
    G0      = config.G0
    FC      = config.FC
    AGM     = config.AGM
    BGM     = config.BGM
    CGM     = config.CGM
    ALPHA1  = config.ALPHA1
    ALPHA2  = config.ALPHA2
    ALPHA3 = config.ALPHA3
    DIF3    = config.DIF3
    DIF4    = config.DIF4
    KP      = config.KP
    C1      = config.C1
    EN1     = config.EN1
    EN3     = config.EN3
    R       = config.R
    WM4     = config.WM4
    WM3     = config.WM3
    WM2     = config.WM2
    WM1     = config.WM1
    HF      = config.HF
    ENMX2   = config.ENMX2
    ENMX3   = config.ENMX3
    TF      = config.TF
    CFL     = config.CFL
    ZEND    = config.ZEND
     
    # Temporary Arrays
    DZ      = 0.0
    Z       = np.zeros(config.LV_MAX)
    H       = np.zeros(config.LV_MAX)
    T       = np.zeros(config.LV_MAX)
    P       = np.zeros(config.LV_MAX)
    DHDZ    = np.zeros(config.LV_MAX)
    DERIV   = np.zeros(config.LV_MAX)

    Z2       = state.Z2
    H2       = state.H2
    T2       = state.T2
    P2       = state.P2
    DHDZ2    = state.DHDZ2
    DERIV2   = state.DERIV2
     
    DERIV   = state.DERIV1
    Z3       = state.Z3
    H3       = state.H3
    T3       = state.T3
    P3       = state.P3
    DHDZ3    = state.DHDZ3
    DERIV3   = state.DERIV3
     
    Z[0]        = Z2[-1]
    H[0]        = H2[-1]
    T[0]        = T2[-1]
    P[0]        = P2[-1]
    DHDZ[0]     = DHDZ2[-1]
    DERIV[0]    = DERIV2[-1]
    GATZ0       = G0 + FC * Z0
    KFLAG       = 0                     # Slope Indicator ( 0: Slope Moving Towards First Peak, 1: Slope has reached First Peak)

     
    TVAP        = UNBAR('TVAP', P[0])
    DELHV       = UNBAR('DHVST', TVAP)
    DELHL       = UNBAR('DHLVST', TVAP)
    
    HL          = (TVAP - TF) * CFL
    HV          = HL + DELHV - DELHL
    
    VP          = UNBAR('TBLVP', T[0])
    CN2H4       = (VP * config.WM4) / (config.R * T[0])
    H4          = UNBAR('TBLH4', T[0])
    AP          = UNBAR('ZTBLAP', Z[0])
    A           = UNBAR('ZTBLA', Z[0])
    
    DZ          = (HL - H[0]) / DHDZ[0] + DZ
    
    LL = 0
    
    C1, C2, C3, C4  = CONC(T[0], P[0], WM4, WM3, WM2, WM1, R, H[0], HF)    
    SUM             = C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4
    FRAC1           = C1 / (WM1*SUM)
    FRAC2           = C2 / (WM2*SUM)
    FRAC3           = C3 / (WM3*SUM)
    FRAC4           = C4 / (WM4*SUM)
    FRAC3D          = (FRAC1/FRAC2 - 1.0) / (3.0 - FRAC1/FRAC2)
    
    print_vapor_data_1(Z[LL], T[LL], P[LL], H[LL], C1, C2, C3, C4) 
    
    CP4   = UNBAR('SHTBL1', T[LL])
    CP3   = UNBAR('SHTBL2', T[LL])
    CP2   = UNBAR('SHTBL3', T[LL])
    CP1   = UNBAR('SHTBL4', T[LL])
    CAV   = (C4*CP4 + C3*CP3 + C2*CP2 + C1*CP1) / (C4 + C3 + C2 + C1)
    WMAV  = (C1+C2+C3+C4) / (C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4)
    H4    = UNBAR('TBLH4', T[LL])
    DELA  = UNBAR('ZTBLD', Z[LL])
    AP    = UNBAR('ZTBLAP', Z[LL])
    A     = UNBAR('ZTBLA', Z[LL])
    
    G, GMMA, K, BETA, DPA = PARAM(T[LL], Z[LL], 1, C4, H4, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, P[LL], KP, C1)
    
    if C4 == 0.0:
        T4 = 0.0
    else:
        DIFN        = DIF4 * ((T[LL]/492.0)**1.823)*14.7/P[LL]
        VIS         = UNBAR('VISVST', T[LL])
        RHO         = P[LL] * WMAV / (R*T[LL])
        AKC         = 0.61 * G / RHO * ((VIS/(RHO*DIFN))**-0.667)*((G/(AP*VIS))**-0.41)
        DERIV3[LL]   = AKC * C4 / DPA
        T4          = AP * DPA * DERIV3[LL]
        
    H3 = UNBAR('TBLH3', T[LL])
    G, GMMA, K, BETA, DPA = PARAM(T[LL], Z[LL], 1, C3, H3, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, P[LL], KP, C1)

    if C3 == 0:
        T3 = 0.0
    else:
        GRAD, TRAD = SGRAD(T[LL], P[LL], G, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3)
        DERIV3[LL] = GRAD/DPA
        T3 = AP * DPA * DERIV3[LL]

    RHOM        = ALPHA3 * C4 * np.exp(-CGM/T[LL])
    T1          = P[LL] * WMAV / (R * T[LL] * G)
    T2          = RHOM * DELA
    DHDZ3[LL]   = -H4/G * (T2+T4) - H3/G * T3 - FC/G * (H[LL]-HF)	
    
    if KFLAG == 1:
        DTDZ    = DHDZ[LL]/CAV
        W1      = C1/RHO
        W2      = C2/RHO
        W3      = C3/RHO
        W4      = C4/RHO
        S1      = 1.0/G
        S5      = FC / (G*RHO)
        DW4DZ   = S1*(FC-T2-T4)-C4*S5
        DW3DZ   = S1*(T2*WM3/WM4+T4*WM3/WM4-T3)-C3*S5
        DW2DZ   = S1*(.5*T2*WM2/WM4+.5*T4*WM2/WM4+.5*T3*WM2/WM3)-C2*S5
        DW1DZ   = S1*(.5*T2*WM1/WM4+.5*T4*WM1/WM4+1.5*T3*WM1/WM3)-C1*S5
        SUMWM   = W1/WM1 + W2/WM2 + W3/WM3 + WM4/WM4
        SMDWDZ  = DW1DZ/WM1 + DW2DZ/WM2 + DW3DZ/WM3 + DW4DZ/WM4
        DMDZ    = -WMAV / SUMWM * SMDWDZ
        DPDZ    = (DELA-1.0) / DELA**3 * (1.75+75.0*VIS*(1.0-DELA)/(A*G))*G**2 / (64.4*A*RHO)
        DPDZ    = DPDZ/144.0
        DRDZR   = DMDZ/WMAV-DTDZ/T[LL] + DPDZ/P[LL]
        T5      = FC/G-DRDZR
        DC4DZ   = T1*(FC-T2-T4)-C4*T5
        DC3DZ   = T1*(T2*WM3/WM4+T4*WM3/WM4-T3)-C3*T5
        DC2DZ   = T1*(.5*T2*WM2/WM4+.5*T4*WM2/WM4+.5*T3*WM2/WM3)-C2*T5
        DC1DZ   = T1*(.5*T2*WM1/WM4+.5*T4*WM1/WM4+1.5*T3*WM1/WM3)-C1*T5
    elif KFLAG == 0:
        Z1      = -H4/(ENMX3 * DHDZ[LL])
        Z2      = 0.05*(ZEND-Z[LL])
        if DHDZ[LL] * (1.0 - Z1/Z2) > 0:
            ...
        else:
            ...
            
        
        
        
        
    
