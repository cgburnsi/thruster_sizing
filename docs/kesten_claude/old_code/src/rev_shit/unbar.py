import numpy as np

def UNBAR(table, data_in):
    """
    Interpolation Routine for Hydrazine and Reaction Products.
    Data is transcribed directly from the FORTRAN 'BLOCK DATA' section
    of NASA G910461-30 (pages 113-115).
    """
    
    table = table.upper()
    calc_val = -1 # Default error value

    # ---
    # Table: Vapor Pressure (psia) vs. Temperature (degR)
    # Source: Page 113, COMMON /FTZ/TBLVP(70) 
    # NOTE: The sample run (Fig. 3c)  uses this table, NOT VPTBL.
    # We must invert the table to get T = f(P).
    # ---
    if table == 'TVAP':
        T_R = [
            492., 519., 528.4, 529.1, 534.6, 534.7, 538.8, 543.9, 545.7, 560.2, 570., 579.3, 
            579.5, 595.3, 610.1, 614.1, 618.1, 627.5, 628.8, 645.7, 650.8, 665.6, 675., 
            686.1, 692.4, 697.5, 744., 798., 852., 942., 1032., 1122., 1176.
        ]
        P_psia = [
            0.052, 0.148, 0.201, 0.207, 0.240, 0.244, 0.282, 0.292, 0.354, 0.545, 0.737, 
            0.973, 0.982, 1.51, 2.20, 2.46, 2.74, 3.41, 3.56, 5.24, 5.97, 8.07, 9.71, 
            11.9, 13.5, 14.7, 33.8, 73.5, 147., 382., 823., 1528., 2131.
        ]
        # Invert the table: Find T as a function of P
        calc_val = np.interp(data_in, P_psia, T_R)  # [degR]

    # ---
    # Table: Enthalpy of Saturated Vapor (Btu/lb) vs. Temperature (degR)
    # Source: Page 115, COMMON /LIZTBL/DHVST(18) [cite: 4582]
    # ---
    elif table == 'DHVST':
        T_R =    [360., 540., 720., 900., 1080., 1260., 1440.]
        DHVST_en = [1390.9, 1332.82, 1280.82, 1279.12, 1237.79, 1208.80, 1139.7]
        calc_val = np.interp(data_in, T_R, DHVST_en) # [Btu/lb]

    # ---
    # Table: Enthalpy of Saturated Liquid (Btu/lb) vs. Temperature (degR)
    # Source: Page 115, COMMON /LIZTBL/DHLVST(18) [cite: 4582]
    # ---
    elif table == 'DHLVST':
        T_R =      [360., 540., 720., 900., 1080.] # Data lines 40 & 60 [cite: 4582]
        DHLVST_en = [100., 360., 534.6, 540.2, 720.] # Data line 60 [cite: 4582]
        calc_val = np.interp(data_in, T_R, DHLVST_en) # [Btu/lb]
        
    # ---
    # Table: Catalyst Particle Radius (ft) vs. Axial Position (ft)
    # Source: Page 22, Fig 2. Sample Data Case 
    # ---
    elif table == 'ZTBLA':
        Z = [0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 
             0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 
             0.1935, 0.2207, 0.2343, 0.2500]
        ZTBLA = [0.001, 0.001, 0.001, 0.001, 0.0064, 0.0064, 0.0064, 0.0064, 
                 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 
                 0.0064, 0.0064, 0.0064, 0.0064, 0.0064]
        calc_val = np.interp(data_in, Z, ZTBLA)  # [ft]

    # ---
    # Table: Particle Surface Area (1/ft) vs. Axial Position (ft)
    # Source: Page 22, Fig 2. Sample Data Case 
    # ---
    elif table == 'ZTBLAP':
        Z = [0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 
             0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 
             0.1935, 0.2207, 0.2343, 0.2500]
        ZTBLAP = [2100.0, 2100.0, 2100.0, 2100.0, 330.0, 330.0, 330.0, 
                  330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 
                  330.0, 330.0, 330.0, 330.0, 330.0]
        calc_val = np.interp(data_in, Z, ZTBLAP)  # [1/ft]

    # ---
    # Table: Void Fraction (-) vs. Axial Position (ft)
    # Source: Page 22, Fig 2. Sample Data Case 
    # ---
    elif table == 'ZTBLD':
        Z = [0.0, 0.0055, 0.0111, 0.0167, 0.0168, 0.0439, 0.0575, 0.0711, 
             0.0847, 0.0983, 0.1119, 0.1255, 0.1391, 0.1527, 0.1663, 0.1799, 
             0.1935, 0.2207, 0.2343, 0.2500]
        ZTBLD = [0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 
                 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34]
        calc_val = np.interp(data_in, Z, ZTBLD)  # [-]

    # ---
    # Table: Vapor Pressure (psia) vs. Temperature (degR)
    # Source: Page 113, COMMON /FTZ/TBLVP(70) 
    # ---
    elif table == 'TBLVP':
        T_R = [
            492., 519., 528.4, 529.1, 534.6, 534.7, 538.8, 543.9, 545.7, 560.2, 570., 579.3, 
            579.5, 595.3, 610.1, 614.1, 618.1, 627.5, 628.8, 645.7, 650.8, 665.6, 675., 
            686.1, 692.4, 697.5, 744., 798., 852., 942., 1032., 1122., 1176.
        ]
        TBLVP_en = [
            0.052, 0.148, 0.201, 0.207, 0.240, 0.244, 0.282, 0.292, 0.354, 0.545, 0.737, 
            0.973, 0.982, 1.51, 2.20, 2.46, 2.74, 3.41, 3.56, 5.24, 5.97, 8.07, 9.71, 
            11.9, 13.5, 14.7, 33.8, 73.5, 147., 382., 823., 1528., 2131.
        ]
        calc_val = np.interp(data_in, T_R, TBLVP_en)  # [psia]

    # ---
    # Table: Heat of Reaction N2H4 (Btu/lb) vs. Temperature (degR)
    # Source: Page 114, COMMON /CCC/H4TBL(40) 
    # ---
    elif table == 'TBLH4':
        T_R = [
            360., 540., 720., 900., 1080., 1260., 1440., 1620., 1800., 1980., 
            2160., 2340., 2520., 2700., 2880., 3060.
        ]
        H4TBL_en = [
            -1951.62, -1919.50, -1896.04, -1895.70, -1882.55, -1878.12, -1879.46,
            -1889.63, -1898.38, -1901.94, -1912.49, -1924.85, -1937.54, -1950.74,
            -1950.74, -1950.74 # Report data list repeats last value
        ]
        calc_val = np.interp(data_in, T_R, H4TBL_en)  # [Btu/lb]

    # ---
    # Table: Heat of Reaction NH3 (Btu/lb) vs. Temperature (degR)
    # Source: Page 114, COMMON /CCC/H3TBL(40) 
    # ---
    elif table == 'H3TBL':
        T_R = [
            360., 540., 720., 900., 1080., 1260., 1440., 1620., 1800., 1980., 
            2160., 2340., 2520., 2700., 2880., 3060.
        ]
        H3TBL_en = [
            1055.57, 1103.97, 1159.35, 1160.40, 1213.46, 1259.64, 1298.00,
            1329.71, 1355.28, 1375.57, 1391.11, 1402.52, 1410.13, 1414.57,
            1416.37, 1416.05
        ]
        calc_val = np.interp(data_in, T_R, H3TBL_en)  # [Btu/lb]

    # ---
    # Table: Specific Heat H2 (Btu/lb-R) vs. Temperature (degR)
    # Source: Page 114, COMMON /DDD/CFTBL1(34) 
    # ---
    elif table == 'CFTBL1':
        T_R = [
            540., 720., 900., 1080., 1260., 1440., 1620., 1800., 
            1980., 2160., 2340., 2520., 2700., 2880., 3060.
        ]
        CFTBL1_en = [
            3.4194, 3.4596, 3.4685, 3.4765, 3.4899, 3.5151, 3.5454, 
            3.5806, 3.6208, 3.6654, 3.7150, 3.7696, 3.8291, 3.8802, 
            3.9288
        ]
        calc_val = np.interp(data_in, T_R, CFTBL1_en)  # [Btu/lb-R]

    # ---
    # Table: Specific Heat N2 (Btu/lb-R) vs. Temperature (degR)
    # Source: Page 114, COMMON /DDD/CFTBL2(34) 
    # ---
    elif table == 'CFTBL2':
        T_R = [
            540., 720., 900., 1080., 1260., 1440., 1620., 1800., 
            1980., 2160., 2340., 2520., 2700., 2880., 3060.
        ]
        CFTBL2_en = [
            0.2485, 0.2495, 0.2524, 0.2569, 0.2624, 0.2682, 0.2738, 
            0.2790, 0.2836, 0.2878, 0.2914, 0.2946, 0.2974, 0.2998, 
            0.3019
        ]
        calc_val = np.interp(data_in, T_R, CFTBL2_en)  # [Btu/lb-R]

    # ---
    # Table: Specific Heat NH3 (Btu/lb-R) vs. Temperature (degR)
    # Source: Page 114, COMMON /DDD/CFTBL3(34) 
    # ---
    elif table == 'CFTBL3':
        T_R = [
            540., 720., 900., 1080., 1260., 1440., 1620., 1800., 
            1980., 2160., 2340., 2520., 2700., 2880., 3060.
        ]
        CFTBL3_en = [
            0.5005, 0.5441, 0.5824, 0.6344, 0.6773, 0.7176, 0.7553, 
            0.7905, 0.8236, 0.8541, 0.8823, 0.9075, 0.9304, 0.9512, 
            0.9697
        ]
        calc_val = np.interp(data_in, T_R, CFTBL3_en)  # [Btu/lb-R]
    
    # ---
    # Table: Specific Heat N2H4 (Btu/lb-R) vs. Temperature (degR)
    # Source: Page 114, COMMON /DDD/CFTBL4(34) 
    # ---
    elif table == 'CFTBL4':
        T_R = [
            540., 720., 900., 1080., 1260., 1440., 1620., 1800., 
            1980., 2160., 2340., 2520., 2700., 2880., 3060.
        ]
        CFTBL4_en = [
            0.5005, 0.5441, 0.5824, 0.6344, 0.6773, 0.7176, 0.7553, 
            0.7905, 0.8236, 0.8541, 0.8823, 0.9075, 0.9304, 0.9512, 
            0.9697
        ]
        calc_val = np.interp(data_in, T_R, CFTBL4_en)  # [Btu/lb-R]

    # ---
    # Tables: Heats of Reaction SHTBL1, SHTBL2, SHTBL4
    # Source: Page 113 
    # ---
    elif table in ['SHTBL1', 'SHTBL2']:
        calc_val = 0.0 # These tables are all zeros
    
    elif table == 'SHTBL4':
        calc_val = -1300.0 # This table is a constant -1300.0
        
    # ---
    # Table: Heat of Reaction NH3 (Btu/lb) vs. Temperature (degR)
    # Source: Page 113, COMMON /FTZ/SHTBL3(34) 
    # ---
    elif table == 'SHTBL3':
        T_R = [
            540., 720., 900., 1080., 1260., 1440., 1620., 1800.,
            1980., 2160., 2340., 2520., 2700., 2880., 3060.
        ]
        SHTBL3_en = [
            -2440., -2455., -2470., -2480., -2488., -2495., -2497., -2495.,
            -2491., -2485., -2474., -2459., -2441., -2419., -2394.
        ]
        calc_val = np.interp(data_in, T_R, SHTBL3_en) # [Btu/lb]
        
    # ---
    # Table: Viscosity (lb/ft-sec) vs. Temperature (degR)
    # Source: Page 113, COMMON /MUVST/VISVST(30) 
    # ---
    elif table == 'VISVST':
        T_R = [
            360.0, 540.0, 720.0, 900.0, 1080.0, 1260.0, 1440.0, 1620.0, 
             1800.0, 1980.0, 2160.0, 2340.0, 2520.0
        ]
        VISVST_en = [
            0.048e-4, 0.070e-4, 0.093e-4, 0.117e-4, 0.141e-4, 0.164e-4, 
            0.186e-4, 0.207e-4, 0.220e-4, 0.247e-4, 0.266e-4, 0.285e-4, 
            0.302e-4
        ]
        calc_val = np.interp(data_in, T_R, VISVST_en)  # [lb/ft-sec]

    else:
        print(f'ERROR - table: {table}, value: {data_in} not calculated correctly')
        calc_val = -1
        
    return calc_val