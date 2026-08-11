from slope02 import slope
import numpy as np

def ZTBLA(x): # [ft] Catalyst Particle Radius
    Z     = [  0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2207, 0.2343, 0.2500]
    ZTBLA = [0.001,  0.001,  0.001,  0.001, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064]
    return np.interp(x, Z, ZTBLA)  

def ZTBLD(x): # [-] Interparticle Void Fraction (DELA)
    Z     = [ 0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2207, 0.2343, 0.2500]
    ZTBLD = [0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34,   0.34]
    return np.interp(x, Z, ZTBLD) 
  
def ZTBLAP(x): # [ft^2] Total Particle External Surface Area
    Z      = [   0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2207, 0.2343, 0.2500]
    ZTBLAP = [2100.0, 2100.0, 2100.0, 2100.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0,  330.0]
    return np.interp(x, Z, ZTBLAP)  

def TBLH4(x): # [?] Not Sure yet what this is
    T     = [     0.0,    180.0,    360.0,    536.4,    540.0,    720.0,    900.0,   1080.0,   1260.0,   1440.0,   1620.0,   1800.0,   1980.0,   2160.0,   2340.0,   2520.0,   2700.0,   2880.0,   3060.0]
    TBLH4 = [-1991.34, -1951.02, -1919.50, -1896.04, -1895.70, -1882.55, -1878.12, -1879.46, -1884.63, -1892.38, -1901.94, -1912.88, -1924.85, -1937.54, -1950.74, -1964.45, -1978.32, -1992.36, -2006.62]
    return np.interp(x, T, TBLH4)

def TBLVP(x): # [psia?] Vapor Pressure of Something? Maybe liquid hydrazine?
    T     = [ 492.0,  519.0, 528.37, 529.08, 534.60, 534.71, 538.84, 543.91, 545.73, 560.20, 569.98, 579.26, 579.48, 595.34, 610.13, 614.08, 618.07, 627.49, 628.82, 645.68, 650.76, 665.57, 674.99 ,686.13, 692.39, 697.47, 744.0, 798.0, 852.0, 942.0, 1032.0, 1122.0, 1176.0]
    TBLVP = [0.0520, 0.1479, 0.2011, 0.2069, 0.2398, 0.2436, 0.2823, 0.2920, 0.3539, 0.5453, 0.7367, 0.9727, 0.9823,  1.510,  2.204,  2.462,  2.740,  3.407,  3.562,  5.240,  5.971,  8.065,  9.711,  11.91,  13.46,  14.70, 33.80, 73.48, 147.0, 382.1,  823.0, 1528.0, 2131.0]
    return np.interp(x, T, TBLVP)

def TVAP(x):
    P_en = np.array([ 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,  550.0,  600.0,  650.0,  700.0,  750.0,  800.0,  850.0,  900.0,  950.0, 1000.0])
    T_en = np.array([770.0, 820.0, 855.0, 880.0, 905.0, 925.0, 945.0, 965.0, 980.0, 995.0, 1010.0, 1025.0, 1035.0, 1050.0, 1060.0, 1070.0, 1080.0, 1090.0, 1100.0, 1110.0])
    return np.interp(x, P_en, T_en)

def DHVST(x):
    TVAP_en  = [ 180.0,   360.0,   534.6,   540.0,   720.0,   900.0,  1080.00]
    DELHV_en = [1390.16, 1332.82, 1280.02, 1279.12, 1237.79, 1208.80, 1189.76]
    return np.interp(x, TVAP_en, DELHV_en)  # [?] Enthalpy of Vaporization 

def DHLVST(x):
    TVAP_en  = [180.0,  360.0,  534.6,  540.0,  720.0,  900.0, 1080.0]
    DELHL_en = [652.14, 665.96, 679.61, 679.89, 700.89, 733.19, 777.22]
    return np.interp(x, TVAP_en, DELHL_en)  # [?] Enthalpy of 

def VISVST(x):
    T_en      = [   360.0,    540.0,    720.0,    900.0,   1080.0,   1260.0, 1440.0, 1620.0, 1800.0, 1980.0, 2160.0, 2340.0, 2520.0]
    VISVST_en = [0.048e-4, 0.070e-4, 0.093e-4, 0.117e-4, 0.141e-4, 0.164e-4, 0.186e-4, 0.207e-4, 0.220e-4, 0.247e-4, 0.266e-4, 0.285e-4, 0.302e-4]
        calc_val = np.interp(data_in, T, VISVST)  # [?] Viscosity Data for N2H4 (maybe N2H4, not sure)


def PARAM(T, ZA, LOP, CC, HR, LVOP, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP):

    if (ZA - Z0) < 0:
        G = G0 + FC * ZA
    else:
        G = GATZ0
        FC = 0.0  # Force FC to 0 in this branch

    if LVOP == 1:
        GMMA = BGM/T
        K = ALPHA2 * np.exp(-GMMA) / C1**1.6
        DPV = DIF3 * (T/492.0)**1.832 * (14.7/PRES) * (1. - np.exp(-.0672*(PRES*492.)/(14.7*T)))
        DPA = DPV
    else:
        GMMA = AGM/T
        K = ALPHA1 * np.exp(-GMMA)
        DPL = DIF4 * (T/492.0)**1.832 * (14.7/PRES) * (1. - np.exp(-.0672*(PRES*492.)/(14.7*T)))
        DPA = DPL

    # The comment questions if it should be CC+HR, not CC*HR. We'll keep CC*HR for accuracy to code.
    BETA = -(CC + HR * DPA) / (KP * T)

    return G, GMMA, K, BETA, DPA

def SLOPE(CG, GMMA, K, BETA, EN12, RATE, DPA, A, DIFF, TEMP=300.0, PRES=1.0, WMAV=1.0, R=0.082, G=9.81, AP=1.0, H=1.0, HL=1.0, C4=1.0, PRINT=False):
    ...


# Initial Concentration of Species in the Vapor Region of the Reactor Bed
def CONC(T, P, WM4, WM3, WM2, WM1, R, H, H4, HF):
    TVAP, PRES = T, P
    XV = -(H - HF) / H4
    C4 = (PRES * WM4) / (R * TVAP) * ((1.0 - XV) / (1.0 + XV))          # Hydrazine Concentration
    C3 = (PRES * WM3) / (R * TVAP) * (XV / (1.0 + XV))                  # Ammonia Concentration
    C2 = (PRES * WM2) / (2.0 * R * TVAP) * (XV / (1.0 + XV))            # Nitrogen Concentration
    C1 = (PRES * WM1) / (2.0 * R * TVAP) * (XV / (1.0 + XV))            # Hydrogen Concentration
    return C1, C2, C3, C4


if __name__ == '__main__':
    
    # Simulation Variables
    max_steps   = 200                       # [-] Maximum Number of Steps Allowed in Calculation
    Z           = np.zeros(max_steps)       # [ft] Axial stations
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
    ALPHA2 = 1.0e10         # [(lb/ft^3)^1.6/sec] Preexponential factor in the rate equation for the catalytic decomposition of ammonia
    EN1 = 1.0               # [-] Order of hydrazine catalytic decomposition reaction with to hydrazine hydrazine

    # Reactor Variables
    G0 = 3.12               # [lb/ft2-s] Inlet Mass Flow Rate
    Z0 = 0.0                # [ft] Axial Distance to the End of a Buried Injector
    FC = 1                  # [lb/ft3-sec] Rate of Feed of Hydrazine into System
    AGM = 2500              # [degR] Activation energy for catalytic decomposition of hydrazine divided by the gas constant
    BGM = 50000             # [degR] Activation energy for catalytic decomposition of ammonia divided by the gas constant
    PRES = 100.0            # [psia] Inlet Chamber Pressure
    KP = 0.4e-4             # [Btu/ft-sec-degR] Thermal Conductivity of the porous catalyst particle (Shell 405)

    # Initialization
    IFC = 1
    TVAP = TVAP(PRES)
    DELHV = DHVST(TVAP)
    DELHL = DHLVST(TVAP)
    HL = (TVAP - TF) * CFL
    HV = HL + DELHV - DELHL
    GATZ0 = G0 + FC * Z0
    H = HF
    
    if FC <= 0: IFC=0
    
    # Main Loop
    i = 1
    while i < max_steps:
        Z[i] = Z[i-1] + DZ
        TEMP = TF + (H - HF) / CFL
        VP = TBLVP(TEMP)
        CN2H4 = (VP * WM4) / (R * TEMP)
        H4 = TBLH4(TEMP)
        AP = ZTBLAP(Z[i])
        A = ZTBLA(Z[i])
                
        G, GMMA, K, BETA, DPA = PARAM(TEMP, Z[i], 1, CN2H4, H4, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP)
        
        RATE, CPA, XOA, MI = slope(CN2H4, GMMA, K, BETA, EN1, DPA, A, DIF4, 
                                   TEMP=300.0, PRES=1.0, WMAV=1.0, R=0.082, G=9.81, AP=1.0, H=1.0, HL=1.0, C4=1.0)
        
        
        # Increment Counter
        i += 1
        
        # Calculate the next DZ
    
    H4 = TBLH4(TVAP)
    print(H4)
    # Post-Processing (Save Data Files, Make Plots, etc.)
    CC1, CC2, CC3, CC4 = CONC(TVAP, PRES, WM4, WM3, WM2, WM1, R, H, H4, HF)
    
    
  