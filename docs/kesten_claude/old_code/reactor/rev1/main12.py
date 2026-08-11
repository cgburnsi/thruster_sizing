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



    
    
def check_relative_convergence(T_prev, T_curr, C_prev, C_curr, tol=0.05):
    """Check if both temperature and concentration are within relative tolerance."""
    temp_diff = abs(T_prev - T_curr) / T_curr
    conc_diff = abs(C_prev - C_curr) / C_curr
    return temp_diff <= tol and conc_diff <= tol

def compute_initial_tps_and_properties(T, H3, KC4, CI4, H4, DCPDX, DP3, HC):
    TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
    if TPS < 0.0:
        TPS = 1.0
    return TPS

def compute_x0_and_rate_params(T, TPS, CI1, CI3, CPS, H3, DP3, KP, BGM, ALPHA2, EN3, A, DCPDX):
    H3P = H3
    DP3P = DP3
    TMTPN = T - TPS
    GAMMA = BGM / TPS
    BETA = -CPS * H3 * DP3 / (KP * TPS)
    K0 = ALPHA2 * np.exp(-GAMMA) * CI1**EN3
    X0 = A - CPS / DCPDX
    X0P = X0
    return H3P, DP3P, TMTPN, GAMMA, BETA, K0, X0, X0P

def clamp_x0_and_recompute_cps(X0, CI3, DP3, A, KC3, T, H4, KC4, CI4, H3, CP3, DCPDX, HC, X0A, NPART, GAMMA, BETA, K0):
    # Clamp X0 and reset profile ratio
    X0 = 0.0
    X0A = 0.0
    
    # Linear profile assumption for CPS
    CPS = CI3 / (DP3 / (A * KC3) + 1.0)
    DCPDX = CI3 / A

    # Recalculate surface temperature TPS
    TPS = T - (H4 * KC4 * CI4 + H3 * CP3 * DCPDX) / HC
    if TPS < 0.0:
        TPS = 1.0

    # Correct CPS using the integral of the profile
    RIESUM = TRAPP(X0A, 1.0, NPART, X0A, CPS, GAMMA, BETA, K0)
    CPS = CI3 - A * RIESUM / KC3

    return X0, X0A, CPS, DCPDX, TPS

def recalculate_cps_with_relaxation(WAF1, WAF2, X0P, X0, KC3, A, DP3, CI3):
    X00 = WAF1 * X0P + WAF2 * X0
    denom = KC3 * A - KC3 * X00
    if abs(denom) < 1e-12:
        raise ZeroDivisionError("Denominator in CPS calculation too small.")
    CPS = CI3 / (1.0 + DP3 / denom)
    return CPS, X00

def adaptive_relaxation_strategy(LP1, WAF1, WAF2, WAF1_CAPPED):
    if LP1 > 25 and not WAF1_CAPPED:
        if WAF1 + 0.05 > 1.0:
            WAF1_CAPPED = True
        else:
            WAF1 += 0.05
            WAF2 = 1.0 - WAF1
    return WAF1, WAF2, WAF1_CAPPED




def EVAL1(A, B):
    return B**3/3.0 - A**3/3.0

def EVAL2(A, B):
    return B**2/2.0 - A**2/2.0
    
def SGRAD(T, P, GG, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3):

    # Initialization
    MAX_ITER            = 25                                # ? Can't recall at the moment
    WAF1                = 0.8                               # [-] Weighting Factor #1
    WAF2                = 1.0 - WAF1                        # [-] Weighting Factor #2
    CPSP                = 0.0                               # Previous estimate for CPS
    TPSPP, TPSP         = 0.0, 0.0                          # Second-to-last estimate and previous estimate for TPS 
    H3P, DP3P           = 0.0, 0.0                          # previous estimate for H3 and DP3
    TMTPN               = 0.0                               # [degR] delta T = gas-surface temperature drop
    CMCPN               = 0.0                               # Concentration Difference between bulk gas and particle surface
    NPART               = 50                                # [-] Number of Steps in the Trapazoidal solver
    LP1                 = 1                                 # [-] Loop Variable
    LP2                 = 1                                 # [-] Loop Variable
    N                   = EN1                               # [?] Set N to EN1 for some reason.  Convient maybe?
    LTFLG               = 0                                 # [-] Flag for something around line 870
    CI1, CI2, CI3, CI4  = C1, C2, C3, C4                    # [lb/ft^3] Interstitial Species Concentrations
    D03, D04            = DIF3, DIF4                        # [ft2/s] Diffusion Coefficients (Gas Phase at STP)
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


    # Locate Suitable X0
    DP3     = DP3F(T, D03, P)
    CPS     = CI3 / 2.0
    CMCPN   = CI3 - CPS
    DCPDX   = KC3/DP3 * (CI3-CPS)
    H3      = UNBAR('H3TBL', T)                             # [BTU/lb] Heat of Reaction for ammonia (Varies with each heat iteration loop)
    H4      = UNBAR('H4TBL', T)                             # [BTU/lb] Heat of Reaction for hydrazine (Constant for each entry to routine)
    
    
    WAF1_CAPPED = False
    converged = False
    
    for LP1 in range(1, 51):
    
        if converged:
            print(f'Converged at LP1 = {LP1}.')
            break
        # (1) - Update TPS, DP3, GAMMA, Betta, K0
        TPS = compute_initial_tps_and_properties(T, H3, KC4, CI4, H4, DCPDX, DP3, HC)
    
        # Update Temperature History
        if LP1 > 1:
            TPSSP = TPSP
            TPSP = TPS
        # (2) Calculate X0
        H3      = UNBAR('H3TBL', TPS)
        DP3     = DP3F(TPS, D03, P)
        H3P, DP3P, TMTPN, GAMMA, BETA, K0, X0, X0P = compute_x0_and_rate_params(T, TPS, CI1, CI3, CPS, H3, DP3, KP, BGM, ALPHA2, EN3, A, DCPDX)

    
        # --- Handle Negative X0 (Linear Profile Fallback) ---
        if X0 < 0.0:
            X0 = 0.0
            X0A = 0.0
            CPS = CI3 / (DP3 / (A * KC3) + 1.0)
            DCPDX = CI3 / A
            print('Negative X0 Calculated during Iteration (Locate Suitable X0 section.)')
    
            TPS = T - (H4 * KC4 * CI4 + H3 * CP3 * DCPDX) / HC
            if TPS < 0.0:
                TPS = 1.0
                print('Negative TPS Calculated During Iteration (Locate Suitable X0 section.)')
    
            # New CPS via TRAPP
            CPSP = CPS
            CMCPO = CMCPN
            print(X0A, 1.0, NPART, X0A, CPS, GAMMA, BETA, K0)
            RIESUM = TRAPP(X0A, 1.0, NPART, X0A, CPS, GAMMA, BETA, K0)
            CPS = CI3 - A * RIESUM / KC3
        else:
            X0A = X0/A
            
    
        # Recalculate CPS using X00
        X00 = WAF1 * X0P + WAF2 * X0
        denom = KC3 * A - KC3 * X00
        if abs(denom) < 1e-12:
            raise ZeroDivisionError("Denominator in CPS calculation is too small")
    
        CPS = CI3 / (1.0 + DP3 / denom)
        CPSP = CPS
        #print(f'[CPS] Recalculated CPS = {CPS:.6f} using X00 = {X00:.6f}')
    
        # --- Adaptive Relaxation ---
        if LP1 > 25 and not WAF1_CAPPED:
            if WAF1 + 0.05 > 1.0:
                print(f'Failed to converge with WAF1={WAF1:.2f}.')
                WAF1_CAPPED = True
            else:
                WAF1 += 0.05
                WAF2 = 1.0 - WAF1
    
        # -----------------------------
        # ✅ Proper Final TP, TMTPN, CMCPN, etc.
        DCPDX = KC3 / DP3 * (CI3 - CPS)
        GRAD = DCPDX * DP3
        AD = HC * (T - TPS)
        TPSPP = TPSP
        P = TPS
        TMTPO = TMTPN
    
        TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
        if TPS <= 0:
            TPS = 1.0
    
        H3 = UNBAR('H3TBL', TPS)
        DP3 = DP3F(TPS, D03, P)
        TMTPN = T - TPS
        GAMMA = BGM / TPS
        BETA = -CPS * H3 * DP3 / (KP * TPS)
        K0 = ALPHA3 * np.exp(-GAMMA) * CI1**EN3
    
        CMCPO = CMCPN
        CMCPN = CI3 - CPS
    
        # ------ Convergence Checks ------
        if check_relative_convergence(TMTPO, TMTPN, CMCPO, CMCPN, tol=0.05):
            if LP1 > 2:
                if TPSP == min(TPS, TPSP, TPSSP) or TPSP == max(TPS, TPSP, TPSSP):
                    TPSP = (TPS + TPSSP) / 2
                    TPSSP = TPSP
                    #print('Smoothed TPSP due to fluctuation in surface temperature')
                else:
                    for LP2 in range(1, 51):

                        # Initialization (Set up spacial grid and integrals)
                        NX      = 24
                        NX1     = NX + 1
                        NXM1    = NX - 1
                        DX      = np.zeros(NX1)
                        CPX     = np.zeros(NX1)
                        CPOX    = np.zeros(NX1)
                        PCPOX   = np.zeros(NX1)
                        RHET    = np.zeros(NX1)
                        
                        VNU     = -KC3 / DP3
                        XA      = X0A                   # Starting Point
                        DELX0A  = (1.0 - X0A) / NX
                        
                        # Linear Approximation (Estimate CP profile with linear shape)
                        for I in range(NX1):
                            if LP2 == 1:
                                CPX[I] = (XA - X0A) / (1.0 - X0A) * CPS  # Linear approximation
                                                
                        # Reaction Rate Evaluation (kinetic model to computer RHET)
                        theta = 1.0 - CPX[I] / CI3
                        denom = 1.0 + BETA * theta
                        exponent = GAMMA * BETA * theta / denom
                        RHET[I] = K0 * CI3**(1 - N) * CPX[I]**N * np.exp(exponent)
                        DX[I] = XA
                        XA += DELX0A
                        
                        # Midpoint Averaging (Smooth CPX, RHET for better integration accuracy)
                        for I in range(NX):
                            CPX[I] = 0.5 * (CPX[I] + CPX[I + 1])
                            RHET[I] = 0.5 * (RHET[I] + RHET[I + 1])
    
                        # Integral Solving: Special case for CPOX[0]
                        DXL     = X0A
                        DXU     = DXL + DELX0A
                        RR1     = 0.0
                        CTRM    = (A * VNU + 1.0) / (A* VNU)
                        
                        for i in range(NX):
                            eval1 = EVAL1(DXL, DXU)
                            eval2 = EVAL2(DXL, DXU)
                            RR1 += RHET[i] * (eval2 - CTRM * eval1)
                            DXL = DXU
                            DXU += DELX0A
                        
                        CPOX[0] = CI3 - A * A / DP3 * RR1
                        if CPOX[0] < 0.0: CPOX[0] = 0.0
                        
                        # Integral Solving (Us EVAL1 and EVAL2) to solve for updated concentration
                        # --- CPOX[1:] General Case Profile ---
                        INT1 = 1
                        K = 1
                        
                        while K < NX:
                            R1 = 0.0
                            R2 = 0.0
                            PS1 = 0.0
                            PS2 = 0.0
                        
                            XA_local = X0A
                            XOA = XA_local + INT1 * DELX0A
                            XA = XOA + DELX0A
                        
                            # First sum (low range integrals)
                            for I in range(INT1):
                                R1 += RHET[I] * EVAL1(XOA, XA)
                                XOA = XA
                                XA += DELX0A
                        
                            R1 *= (1.0 / XOA - CTRM)
                        
                            # Second sum (upper range integrals)
                            XAD = XA
                            XA -= DELX0A
                            for I in range(INT1, NXM1):
                                PS1 += RHET[I+1] * EVAL2(XA, XAD)
                                PS2 += RHET[I+1] * EVAL1(XA, XAD)
                                XA = XAD
                                XAD += DELX0A
                        
                            R2 = PS1 - CTRM * PS2
                            CPOX[K] = CI3 - (A**2 / DP3) * (R1 + R2)
                            if CPOX[K] < 0.0:
                                CPOX[K] = 0.0
                        
                            INT1 += 1
                            K += 1
                            
                        # Boundary Profile Calculation (Handle X=0 and X=A specifically)
                        DXL = X0A
                        DXU = DXL + DELX0A
                        RR2 = 0.0
                        
                        for I in range(NX):
                            RR2 += RHET[I] * EVAL1(DXL, DXU)
                            DXL = DXU
                            DXU += DELX0A
                        
                        CPOX[NX] = CI3 - (A**2 / DP3) * (1.0 - CTRM) * RR2
                        if CPOX[NX] < 0.0:
                            CPOX[NX] = 0.0
    
                        # Update Temperature Related Quantities Before Checking Convergence
                        # --- Update DCPDX, TPS, and DP3 ---
                        DCPDX = KC3 / DP3 * (CI3 - CPOX[NX])
                        
                        H3P   = H3
                        DP3P  = DP3
                        
                        TPS   = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
                        if TPS < 0.0:
                            TPS = 1.0  # Clamp for stability
                        
                        H3    = UNBAR('H3TBL', TPS)
                        DP3   = DP3F(TPS, D03, P)
                        
                        TMTPO = TMTPN
                        TMTPN = T - TPS

                        # Convergence Check (Temperature and Concentration Check for Stopping)
                        if LP2 > 1:
                            CMCPO = CMCPN
                            CMCPN = CI3 - CPOX[NX]
                            
                            if check_relative_convergence(TMTPN, TMTPO, CMCPN, CMCPO, tol=0.05):
                                print(f'Converged at LP2 = {LP2}.')
                                break
                        
                        # Profile Blending (Smooth new and old profiles to stabilize convergence)
                        for I in range(NX1):
                            if LP2 % 5 == 0:
                                CPX[I] = 0.5 * (CPOX[I] + PCPOX[I])
                            else:
                                CPX[I] = 0.8 * CPX[I] + 0.2 * CPOX[I]
                            
                            PCPOX[I] = CPOX[I]
                        
                        # Final Gradient Calculation (Derive DCDPX, GRAD, and TGRAD)
                        CMCPN = CI3 - CPX[NX1 - 1]
                        DCPDX = KC3 / DP3P * (CI3 - CPX[NX1 - 1])
                        
                        TPS = T - (H4 * KC4 * CI4 + H3P * DP3P * DCPDX) / HC
                        if TPS < 0.0:
                            TPS = 1.0
                        
                        H3 = UNBAR("H3TBL", TPS)
                        DP3 = DP3F(TPS, D03, P)
                        
                        TMTPO = TMTPN
                        TMTPN = T - TPS
                        
                        GRAD = DCPDX * DP3
                        TGRAD = HC * (T - TPS)

                        #print(f"KC3 * (CI3 - CPS) = {GRAD:.6e}")
                        #print(f"HC  * (T  - TPS) = {TGRAD:.6e}")
                        
                        # Update temperature and concentration histories for next iteration
                        TMTPO = TMTPN
                        CMCPO = CMCPN
                        
                    else:
                        # Only Execute if LP2 loop runs 50 times without breaking
                        #print(f"Unable to converge on CPS in 50 tries... CP(X/A) = {CPOX[NX1 - 1]:.6e}")
                        #print(f"KC3*(CI3 - CPS) = {GRAD:.6e}\nHC*(T - TPS) = {TGRAD:.6e}")
                        break  # Exit LP1 loop as well
                        
                    converged = True
    else:
        LTFLG = 1
        print('❌ Failed to converge after 50 iterations.')
    
    return GRAD, TGRAD



def REDIVD(DZ1, DTDZ, NINT, JFLAG, I, LL, Z, Z0, TEMP):
    """
    Translates Fortran subroutine REDIVD to Python.

    Inputs:
        DZ1   - Initial DZ step
        DTDZ  - dT/dZ at current LL
        NINT  - Current interval counter
        JFLAG - DZ condition flag (to be updated)
        I     - Nesting index (to be updated)
        LL    - Current axial index
        Z     - Array of axial positions
        Z0    - Injector start position
        TEMP  - Current temperature at LL

    Outputs:
        Updated values of (NINT, JFLAG, I, DZ)
    """
    I = 0
    NESTCT = 1

    if NINT > 1:
        NESTCT = NINT

    while True:
        I += 1
        XSIZE = 2 ** I
        NINT = int(XSIZE)
        DZ = DZ1 / XSIZE

        # Condition 1: If predicted temp change is acceptable
        if abs(DTDZ) * DZ > 0.01 * TEMP:
            continue  # keep dividing

        # Condition 2: If DZ is safe with respect to Z and Z0
        delta_check = 1.0 + DZ / (Z[LL] - Z0) + 0.01 * Z0 / abs(Z[LL] - Z0)
        if delta_check < 0.0:
            continue  # keep dividing

        # Both conditions are satisfied
        break

    NINT *= NESTCT
    JFLAG = 1

    return DZ, NINT, JFLAG, I





if __name__ == '__main__':
    
    # Initialization
    config  = SimConfig()
    state   = SimState()
    
    print('----- LIQUID ZONE -------------')
    state.Z1, state.T1, state.P1, state.H1, state.DHDZ1, state.DERIV1 = LIQUID(state, config)    
    print('----- LIQUID - VAPOR ZONE -----')
    state.Z2, state.T2, state.P2, state.H2, state.DHDZ2, state.DERIV2 = LIQUIDVAPOR(state, config)
    print('----- VAPOR ZONE --------------')
    
    
    
    # Step 1: Start of Vapor Phase Calculations 
    
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
    #H       = np.zeros(config.LV_MAX)
    #T       = np.zeros(config.LV_MAX)
    #P       = np.zeros(config.LV_MAX)
    DHDZ    = np.zeros(config.LV_MAX)
    #DERIV   = np.zeros(config.LV_MAX)

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
    
    Z[0]     = Z2[-1]
    H        = H2[-1]
    T        = T2[-1]
    P        = P2[-1]
    DHDZ[0]  = DHDZ2[-1]
    DERIV    = DERIV2[-1]
    GATZ0       = G0 + FC * Z0
    KFLAG       = 0                     # Slope Indicator ( 0: Slope Moving Towards First Peak, 1: Slope has reached First Peak)

    DTDZ        = 0
    
    TVAP        = UNBAR('TVAP', P)
    DELHV       = UNBAR('DHVST', TVAP)
    DELHL       = UNBAR('DHLVST', TVAP)
    
    TEMP        = TVAP
    PRES        = P
    
    HL          = (TVAP - TF) * CFL
    HV          = HL + DELHV - DELHL
    
    VP          = UNBAR('TBLVP', TEMP)
    CN2H4       = (VP * config.WM4) / (config.R * TEMP)
    H4          = UNBAR('TBLH4', TEMP)
    AP          = UNBAR('ZTBLAP', Z[0])
    A           = UNBAR('ZTBLA', Z[0])
    
    ZBOUND      = ZEND
    DZ          = (HL - H) / DHDZ[0] + DZ
    
    # Variables used in the Peak temperature checks below
    KOUNT   = 1
    N       = 5
    THREE   = 3.0
    JFLAG   = 0
    IFC     = 0
    NINT    = 0
    I       = 0
    PFRC3D = 0.0  # Previous FRAC3D, used to check slope changes
    
    
    # Step 2: Initialize Parameters for Vapor Region
    C1, C2, C3, C4  = CONC(TEMP, PRES, WM4, WM3, WM2, WM1, R, H, HF)    
    SUM             = C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4
    FRAC1           = C1 / (WM1*SUM)
    FRAC2           = C2 / (WM2*SUM)
    FRAC3           = C3 / (WM3*SUM)
    FRAC4           = C4 / (WM4*SUM)
    FRAC3D          = (FRAC1/FRAC2 - 1.0) / (3.0 - FRAC1/FRAC2)
    
    LL = 0
    while Z[LL] < ZEND and LL + 1 < config.LV_MAX:
        
        # Output to Console the Current Data
        print_vapor_data_1(Z[LL], TEMP, PRES, H, C1, C2, C3, C4) 

        
        # Step 3: Fluid Properties Lookup and/or Calculations
        CP4   = UNBAR('SHTBL1', TEMP)
        CP3   = UNBAR('SHTBL2', TEMP)
        CP2   = UNBAR('SHTBL3', TEMP)
        CP1   = UNBAR('SHTBL4', TEMP)
        CAV   = (C4*CP4 + C3*CP3 + C2*CP2 + C1*CP1) / (C4 + C3 + C2 + C1)
        WMAV  = (C1+C2+C3+C4) / (C1/WM1 + C2/WM2 + C3/WM3 + C4/WM4)
        H4    = UNBAR('TBLH4', TEMP)
        DELA  = UNBAR('ZTBLD', Z[LL])
        AP    = UNBAR('ZTBLAP', Z[LL])
        A     = UNBAR('ZTBLA', Z[LL])
        
        G, GMMA, K, BETA, DPA = PARAM(TEMP, Z[LL], 1, C4, H4, 0, Z0, G0, FC, GATZ0, AGM, BGM, 
                                      ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1)
        
        
        # Step 4:   Calculate GRAD and TGRAD.  There are several parts to this step
        #           1) Check if hydrazine is present.  
        #           2) Check if ammonia is present
        
        # Step 4-1: Check if hydrazine is present.  If C4 = 0, then no hydrazine is at the current axial station
        T4 = 0.0        # [lb/ft^2-s] Mass Flux per unit projected area?  Maybe that's correct??
        if C4 != 0.0:
            DIFN        = DIF4 * ((TEMP/492.0)**1.823)*14.7/PRES
            VIS         = UNBAR('VISVST', TEMP)
            RHO         = PRES * WMAV / (R*TEMP)
            AKC         = 0.61 * G / RHO * ((VIS/(RHO*DIFN))**-0.667)*((G/(AP*VIS))**-0.41)
            DERIV       = AKC * C4 / DPA
            T4          = AP * DPA * DERIV
        
        H3 = UNBAR('TBLH3', TEMP)
        G, GMMA, K, BETA, DPA = PARAM(TEMP, Z[LL], 1, C3, H3, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1)

        # Step 4-2: Check if ammonia is present.  If C3 = 0, then no ammonia is at the current axial station and skips SGRAD completely
        T3 = 0.0
        if C3 != 0.0:
            GRAD, TGRAD = SGRAD(TEMP, PRES, G, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3)
            DERIV       = GRAD/DPA
            T3          = AP * DPA * DERIV
        
        # 4: Main Gradient Calculations. Compute DHDZ
        RHOM        = ALPHA3 * C4 * np.exp(-CGM/TEMP)
        T1          = PRES * WMAV / (R * TEMP * G)
        T2          = RHOM * DELA
        DHDZ[LL]    = -H4/G * (T2+T4) - H3/G * T3 - FC/G * (H-HF)	
        
        # Step 5:   Determine if we are at the Peak Temperature. There are several parts to this step        
        if KFLAG == 0:    # Slope is moving to the peak
            Z1      = -H4/(ENMX3 * DHDZ[LL])
            Z2      = 0.05*(ZEND-Z[LL])
            
            if DHDZ[LL] * (1.0 - Z1/Z2) > 0:
                DZ = Z1
            else:
                print('nonono')
                DELTAZ = ZEND - Z[LL]
                ZBOUND = ZEND + DELTAZ / 3.0
                KFLAG = 1
                KOUNT += 1
        else:               # KFLAG == 1 - Slope as reached the peak
            pass        # Fall through to the next part of the code
        
        # Step 6:   Compute DTDZ and other derivatives
        DTDZ = DHDZ[LL] / CAV
        
        # Species mass fractions
        W1      = C1/RHO
        W2      = C2/RHO
        W3      = C3/RHO
        W4      = C4/RHO
        
        # Intermediate factors
        S1      = 1.0/G
        S5      = FC / (G*RHO)
            
        # Species mass fraction gradients (DWnDZ)
        DW4DZ   = S1*(FC-T2-T4)-C4*S5
        DW3DZ   = S1*(T2*WM3/WM4+T4*WM3/WM4-T3)-C3*S5
        DW2DZ   = S1*(.5*T2*WM2/WM4+.5*T4*WM2/WM4+.5*T3*WM2/WM3)-C2*S5
        DW1DZ   = S1*(.5*T2*WM1/WM4+.5*T4*WM1/WM4+1.5*T3*WM1/WM3)-C1*S5
        
        # Summation terms for average molecular weight and its gradient
        SUMWM   = W1/WM1 + W2/WM2 + W3/WM3 + WM4/WM4
        SMDWDZ  = DW1DZ/WM1 + DW2DZ/WM2 + DW3DZ/WM3 + DW4DZ/WM4
            
        # === Pressure gradient (Fortran lines 930–950) ===
        DPDZ = ((DELA - 1.0) / DELA**3) * (1.75 + 75.0 * VIS * (1.0 - DELA) / (A * G)) * G**2
        DPDZ /= (64.4 * A * RHO)
        DPDZ /= 144.0  # convert to appropriate units
        
        # Gradient of average molecular weight
        DMDZ    = -WMAV / SUMWM * SMDWDZ
        
        # === Density gradient and auxiliary terms (lines 960–970) ===
        DRDZR = DMDZ / WMAV - DTDZ / TEMP + DPDZ / PRES
        T5 = FC / G - DRDZR
        
        # === Species concentration gradients (lines 980–1010) ===
        DC4DZ = T1 * (FC - T2 - T4) - C4 * T5
        DC3DZ = T1 * (T2 * WM3 / WM4 + T4 * WM3 / WM4 - T3) - C3 * T5
        DC2DZ = T1 * (0.5 * T2 * WM2 / WM4 + 0.5 * T4 * WM2 / WM4 + 0.5 * T3 * WM2 / WM3) - C2 * T5
        DC1DZ = T1 * (0.5 * T2 * WM1 / WM4 + 0.5 * T4 * WM1 / WM4 + 1.5 * T3 * WM1 / WM3) - C1 * T5

        # Step 7: Re-check if TEMP vs. Z has reached the initial peak
        #H += DHDZ[LL] * DZ
        #if H < HV:
        #    print('Liquid Hydrazine at vapor interface')
        #    break
        
        # If still in pre-peak region
        if KFLAG == 0:
            pass  # Already updated H and checked for HV
        
        # Post-slope-peak logic
        else:
            if JFLAG == 1:  # DZ too small initially
                if abs(DTDZ) * DZ > 0.01 * TEMP:
                    DZ1 = DZ
                    DZ, NINT, JFLAG, I = REDIVD(DZ1, DTDZ, NINT, JFLAG, I, LL, Z, Z0, TEMP)
                    print(f"\n KOUNT={KOUNT:2d} --- THIS INTERVAL HAS BEEN REDIVIDED {NINT:4d} TIMES")
                else:
                    # Check if we're near the end of the injector
                    injector_end = (1.0 + DZ / (Z[LL] - Z0) + 0.01 * Z0 / abs(Z[LL] - Z0)) > 0.0
                    if not injector_end:
                        DZ1 = DZ
                        DZ, NINT, JFLAG, I = REDIVD(DZ1, DTDZ, NINT, JFLAG, I, LL, Z, Z0, TEMP)
                        print(f"\n KOUNT={KOUNT:2d} --- THIS INTERVAL HAS BEEN REDIVIDED {NINT:4d} TIMES")
        
            else:  # JFLAG == 0: DZ increment is acceptable
                if KOUNT in (4, 6, 8, 10, 12, 14):
                    N -= 1
                KOUNT += 1
                DZ = DELTAZ / (THREE ** N)
        
                if FC <= 0.0 and IFC == 0:
                    pass  # H already updated above, no further checks
                elif abs(DTDZ) * DZ > 0.01 * TEMP:
                    DZ1 = DZ
                    DZ, NINT, JFLAG, I = REDIVD(DZ1, DTDZ, NINT, JFLAG, I, LL, Z, Z0, TEMP)
                    print(f"\n KOUNT={KOUNT:2d} --- THIS INTERVAL HAS BEEN REDIVIDED {NINT:4d} TIMES")
                else:
                    # Check again if we are near the end of the injector
                    injector_end = (1.0 + DZ / (Z[LL] - Z0) + 0.01 * Z0 / abs(Z[LL] - Z0)) > 0.0
                    if not injector_end:
                        DZ1 = DZ
                        DZ, NINT, JFLAG, I = REDIVD(DZ1, DTDZ, NINT, JFLAG, I, LL, Z, Z0, TEMP)
                        print(f"\n KOUNT={KOUNT:2d} --- THIS INTERVAL HAS BEEN REDIVIDED {NINT:4d} TIMES")
        
        # Final H check after DZ refinement
        H += DHDZ[LL] * DZ
        if H < HV:
            print('Liquid Hydrazine - final recheck')
            break
               
        # Step 12: Perform Calculations for Temp, Press, Cx, and mole fractions for Cx
        TEMP += DTDZ * DZ
        PRES += DPDZ * DZ
        C4   += DC4DZ * DZ
        C3   += DC3DZ * DZ
        C2   += DC2DZ * DZ
        C1   += DC1DZ * DZ
        
        SUM = C1 / WM1 + C2 / WM2 + C3 / WM3 + C4 / WM4
        if C4 < 0.0:
            SUM -= C4 / WM4
        if C3 < 0.0:
            SUM -= C3 / WM3
        
        FRAC1 = C1 / (WM1 * SUM)
        FRAC2 = C2 / (WM2 * SUM)
        FRAC3 = C3 / (WM3 * SUM)
        FRAC4 = C4 / (WM4 * SUM)
        
        FRAC3D = (FRAC1 / FRAC2 - 1.0) / (3.0 - FRAC1 / FRAC2)

        # Step 13: Determine if FRAC3D is too large
        if KFLAG == 1:
            PFRC3D = FRAC3D
            
            if C4 < 0.0:
                C4 = 0.0
                FRAC4 = 0.0
            
            if C3 < 0.0:
                C3 = 0.0
                FRAC3 = 0.0
        else:               # KFLAG == 0
            # If the relative difference of successive FRAC3D's is greater than 5 percent,
            # we recalculate with smaller DZ increment
            if (FRAC3D-PFRC3D) < 0.05:
                PFRC3D = FRAC3D
                
                if C4 < 0.0:
                    C4 = 0.0
                    FRAC4 = 0.0
                
                if C3 < 0.0:
                    C3 = 0.0
                    FRAC3 = 0.0
            else:
                H       = H-DHDZ[LL]*DZ
                TEMP    = TEMP - DPDZ*DZ
                PRES    = PRES - DPDZ*DZ
                C4      = C4 - DC4DZ*DZ
                C3      = C3 - DC3DZ*DZ
                C2      = C2 - DC2DZ*DZ
                C1      = C1 - DC1DZ*DZ
                DZ=DZ/2.0
                
                H = H + DHDZ[LL]*DZ
                if H < HV:
                    print('Liquid Hydrazine at the liquid-vapor to vapor interface')
                    break
                else:
                    pass
       
        LL += 1
        Z[LL] = Z[LL-1] + DZ
        if JFLAG == 1: NINT = NINT=1
        if NINT == 0: JFLAG = 0
        
        if Z[LL] > ZBOUND:
            if IFC == 1:
                print('Z from Vapor Region are printed here')
                MBSS = (C1+C2+C3+C4)/(C1/WM1+C2/WM2+C3/WM3+C4/WM4)
                print('steady state values for mbar and g at the end of the bed (MBSS, G)')
            else:
                pass
                
        if KOUNT > 15 and JFLAG == 0:
            if IFC == 1:
                print('Z from Vapor Region are printed here')
                MBSS = (C1+C2+C3+C4)/(C1/WM1+C2/WM2+C3/WM3+C4/WM4)
                print('steady state values for mbar and g at the end of the bed (MBSS, G)')
            else:
                pass
        else:
            INJECT=0
            
        
        
        
    
