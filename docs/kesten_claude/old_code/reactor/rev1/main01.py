import numpy as np
from slope01 import slope


def TVAP_LUT(P_input):
    P_en = np.array([ 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,  550.0,  600.0,  650.0,  700.0,  750.0,  800.0,  850.0,  900.0,  950.0, 1000.0])
    T_en = np.array([770.0, 820.0, 855.0, 880.0, 905.0, 925.0, 945.0, 965.0, 980.0, 995.0, 1010.0, 1025.0, 1035.0, 1050.0, 1060.0, 1070.0, 1080.0, 1090.0, 1100.0, 1110.0])
    return np.interp(P_input, P_en, T_en)

def DELHV_LUT(TVAP_input):
    TVAP_en  = [ 180.0,   360.0,   534.6,   540.0,   720.0,   900.0,  1080.00]
    DELHV_en = [1390.16, 1332.82, 1280.02, 1279.12, 1237.79, 1208.80, 1189.76]
    return np.interp(TVAP_input, TVAP_en, DELHV_en)  # [?] Enthalpy of Vaporization 

def DELHL_LUT(TVAP_input):
    TVAP_en  = [180.0,  360.0,  534.6,  540.0,  720.0,  900.0, 1080.0]
    DELHL_en = [652.14, 665.96, 679.61, 679.89, 700.89, 733.19, 777.22]
    return np.interp(TVAP_input, TVAP_en, DELHL_en)  # [?] Enthalpy of 

def VP_LUT(TEMP_input):
    T_en     = np.array([492.0, 519.0, 528.37, 529.08, 534.60, 534.71, 538.84, 543.91, 545.73, 560.20, 569.98, 579.26, 579.48, 595.34, 610.13, 614.08, 618.07, 627.49, 628.82, 645.68, 650.76, 665.57, 674.99 ,686.13, 692.39, 697.47, 744.0, 798.0, 852.0, 942.0, 1032.0, 1122.0, 1176.0])
    TBLVP_en = np.array([0.0520, 0.1479, 0.2011, 0.2069, 0.2398, 0.2436, 0.2823, 0.2920, 0.3539, 0.5453, 0.7367, 0.9727, 0.9823, 1.510, 2.204, 2.462, 2.740, 3.407, 3.562, 5.240, 5.971, 8.065, 9.711, 11.91, 13.46, 14.70, 33.80, 73.48, 147.0, 382.1, 823.0, 1528.0, 2131.0])
    return np.interp(TEMP_input, T_en, TBLVP_en)

def AP_LUT(Z_input):
    Z_en      = np.array([0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2207, 0.2343, 0.2500])
    ZTBLAP_en = np.array([2100.0, 2100.0, 2100.0, 2100.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0])
    return np.interp(Z_input, Z_en, ZTBLAP_en)
    
def A_LUT(Z_input):
    Z_en     = np.array([0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2207, 0.2343, 0.2500])
    ZTBLA_en = np.array([0.001, 0.001, 0.001, 0.001, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064])
    return np.interp(Z_input, Z_en, ZTBLA_en)

def H4_LUT(TEMP_input):
    T_en     = np.array([0.0, 180.0, 360.0, 536.4, 540.0, 720.0, 900.0, 1080.0, 1260.0, 1440.0, 1620.0, 1800.0, 1980.0, 2160.0, 2340.0, 2520.0, 2700.0, 2880.0, 3060.0])
    TBLH4_en = np.array( [-1991.34, -1951.02, -1919.50, -1896.04, -1895.70, -1882.55, -1878.12, -1879.46, -1884.63, -1892.38, -1901.94, -1912.88, -1924.85, -1937.54, -1950.74, -1964.45, -1978.32, -1992.36, -2006.62])
    return np.interp(TEMP_input, T_en, TBLH4_en)  # [?] 






    
# Inputs
nofz = 20               # [-] Number of Axial Stations (Z's) to be used in the tables
Z0 = 0                  # [ft] Axial Distance to the End of a Buried Injector
zend = 0.25             # [ft] Catalyst Bed Lengt
G0 = 3.12               # [lb/ft2-s] Inlet Mass Flow Rate
TF = 530                # [degR] Temperature of liquid hydrazine entering the bed
PRES = 100.0            # [psia] Inlet Chamber Pressure
FC = 1                  # [lb/ft3-sec] Rate of Feed of Hydrazine into System
HF = 0                  # [Btu/lb] Enthalpy of liquid hydrazine entering the bed
ALPHA1 = 1.0e10         # [1/sec] Preexponetial factor in the rate equation for the catalytic decomposition of hydrazine
alpha2 = 1.0e10         # [(lb/ft^3)^1.6/sec] Preexponential factor in the rate equation for the catalytic decomposition of ammonia
alpha3 = 2.14e10        # [1/sec] Preexponetial factor in the rate equation for the thermal decomposition of hydrazine
EN1 = 1.0               # [-] Order of hydrazine catalytic decomposition reaction with to hydrazine hydrazine
en2 = 1.0               # [-] Order of ammonia catalytic decomposition reaction with respect to ammonia
en3 = -1.6              # [-] Order of ammonia catalytic decomposition reaction with respect to hydrogen
AGM = 2500              # [degR] Activation energy for catalytic decomposition of hydrazine divided by the gas constant
BGM = 50000             # [degR] Activation energy for catalytic decomposition of ammonia divided by the gas constant
cgm = 33000             # [degR] Activation energy for thermal decomposition of hydrazine divided by the gas constant
R = 10.73               # [psia-ft3/lb-mol-degR] Gas Constant
dif3 = 0.17e-3          # [ft2/s] Diffusion coefficient of ammonia in the gas phase at STP
DIF4 = 0.95e-4          # [ft2/s] Diffusion coefficient of hydrazine in the gas phase at STP
wm1 = 2.016             # [lb/lb-mol] Molecular Weight of Hydrogen
wm2 = 28.016            # [lb/lb-mol] Molecular Weight of Nitrogen
wm3 = 17.032            # [lb/lb-mol] Molecular Weight of Amonia
WM4 = 32.048            # [lb/lb-mol] Molecular Weight of Hydrazine
KP = 0.4e-4             # [Btu/ft-sec-degR] Thermal Conductivity of the porous catalyst particle (Shell 405)
CFL = 0.7332            # [Btu/lb-degR] Specific Heat of Liquid Hydrazine
ENMX1 = 200             # [-] Constant used to determine axial station increments in liquid region
enmx2 = 40              # [-] constant used to determine axial station increments in the liquid-vapor region
enmx3 = 80              # [-] constant used to determine axial station increments in the vapor region

TBLVP = np.zeros(70)
TBLH4 = np.zeros(42)
TBLH3 = np.zeros(42)
SHTBL1 = np.zeros(34)
SHTBL2 = np.zeros(34)
SHTBL3 = np.zeros(34)
SHTBL4 = np.zeros(34)
ZTBLD = np.zeros(46)
ZTBLAP = np.zeros(46)
ZTBLA = np.zeros(46)
VISVST = np.zeros(30)
DHVST = np.zeros(18)
DHLVST = np.zeros(18)
VPTBL = np.zeros(44)
IFC = 1                         # main 220

DERIV = np.zeros(250)
DHDZ = np.zeros(250)
Z = np.zeros(25)


nztbl = 2*nofz+4
nofz4 = nofz+4
nofz5 = nofz4+1


if __name__ == '__main__':

    VPTBL[0]  = TVAP  = TVAP_LUT(PRES)
    DHVST[0]  = DELHV = DELHV_LUT(820)
    DHLVST[0] = DELHL = DELHL_LUT(820)
    HL        = (TVAP - TF) * CFL
    HV        = HL + DELHV - DELHL
    GATZ0     = G0 + FC * Z0
    if GATZ0 == 0: 
        IFC = 0

    #Z_station = np.array([0.0, 0.55e-2, 0.111e-1, 0.167e-1, 0.168e-1, 0.439e-1, 0.575e-1, 0.711e-1, 0.847e-1, 0.983e-1, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 0.1935, 0.2071, 0.2207, 0.2207, 0.2343, 0.25])
    #ZTBLA     = np.array([0.1e-2, 0.1e-2, 0.1e-2, 0.1e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2, 0.64e-2])
    #ZTBLAP    = np.array([.21e4, .21e4, 0.21e4, 0.21e4, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3, 0.33e3])
    #ZTBLD     = np.array([0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0, 0.34e0])


    print('--- Entering Liquid Region ---')
    MFLAG = 0
    DZ = 0.0
    Z[0] = 0
    H = HF
    II = 1
    
    Z[II] = Z[II-1] + DZ
    TEMP = TF + (H - HF) / CFL
    VP = VP_LUT(TEMP)
    CN2H4 = (VP * WM4) / (R*TEMP)
    H4 = H4_LUT(TEMP)
    AP = AP_LUT(Z[II])
    A = A_LUT(Z[II])
    
    # PARAM SUBROUTINE (   T,    ZA, LOP,    CC, HR, LVOP, G, GMMA, K, DPA, BETA)
    # PARAM SUBROUTINE (TEMP, Z(II),   1, CN2H4, H4,    0, G, GMMA, K, DPA, BETA)
    # This section of code was the param subroutine.  I've put it inline here.
    if (Z[II]-Z0) < 0.0:
        G = G0 + FC*Z[II]
    else:
        G = GATZ0
        FC = 0.0
        
    GMMA = AGM/TEMP
    # Calculate K, DPA for N2H4
    K = ALPHA1*np.exp(-GMMA)
    DPL=DIF4*(TEMP/492.0)**1.832 * (14.7/PRES)*(1.0-np.exp(-0.0672*(PRES*492.0)/(14.7*TEMP)))
    DPA=DPL
    BETA=-(CN2H4*H4*DPA)/(KP*TEMP)
    # End of PARAM subroutine

    RATE, CPA, XOA, MI = slope(CN2H4, GMMA, K, BETA, EN1, DERIV[II], DPA, A, DIF4, TEMP=300.0, PRES=1.0, 
              WMAV=1.0, R=0.082, G=9.81, AP=1.0, H=1.0, HL=1.0, C4=1.0, PRINT=False)
    
    if H-HL == 0:
        if MI>20: DERIV[II] = DERIV[II-1]
    else:
        DHDZ[II] = -(H4*DPA*AP*DERIV[II]+FC*(H-HF))/G
        DZ = -H4/(ENMX1*DHDZ[II])
        print(f'Z = {Z[II]}, TEMP = {TEMP}, H = {H}, DHDZ = {DHDZ}')
        
    if (H-HL) != 0:
        H = H+DHDZ[II]*DZ
        
    if (H-HL) < 0:
        II = II+1
    
        