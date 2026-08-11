from print_data import print_lqvp_data
import numpy as np
from unbar import UNBAR
from param2 import PARAM


def LIQUIDVAPOR(state, config):
    # Retrieve Config Parameters
    Z0      = config.Z0
    G0      = config.G0
    FC      = config.FC
    AGM     = config.AGM
    BGM     = config.BGM
    ALPHA1  = config.ALPHA1
    ALPHA2  = config.ALPHA2
    DIF3    = config.DIF3
    DIF4    = config.DIF4
    KP      = config.KP
    C1      = config.C1
    EN1     = config.EN1
    R       = config.R
    WM4     = config.WM4
    HF      = config.HF
    ENMX2   = config.ENMX2
    TF      = config.TF
    CFL     = config.CFL
    
    # Temporary Arrays
    Z       = np.zeros(config.LV_MAX)
    H       = np.zeros(config.LV_MAX)
    T       = np.zeros(config.LV_MAX)
    P       = np.zeros(config.LV_MAX)
    DHDZ    = np.zeros(config.LV_MAX)
    DERIV   = np.zeros(config.LV_MAX)
    
    DERIV   = state.DERIV1
    Z2       = state.Z1
    H2       = state.H1
    T2       = state.T1
    P2       = state.P1
    DHDZ2    = state.DHDZ1
    DERIV2   = state.DERIV1
    
    DZ          = 0.0
    Z[0]        = Z2[-1]
    H[0]        = H2[-1]
    T[0]        = T2[-1]
    P[0]        = P2[-1]
    DHDZ[0]     = DHDZ2[-1]
    DERIV[0]    = DERIV2[-1]
    
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

    JJ = 1
    while JJ < config.LV_MAX:
        Z[JJ] = Z[JJ - 1] + DZ
        T[JJ] = T[JJ-1]
    
        P[JJ] = P[JJ - 1]

        TEMP = TVAP
        DERIV[JJ] = DERIV[JJ - 1]

        # Use previous Z for lookups
        Z_prev = Z[JJ - 1]
        H4 = UNBAR('TBLH4', TEMP)
        AP = UNBAR('ZTBLAP', Z_prev)
        A  = UNBAR('ZTBLA', Z_prev)

        GATZ0 = G0 + FC * Z0
        G, GMMA, K, BETA, DPA = PARAM(TEMP, Z_prev, 1, 0.0, 0.0, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, P[0], KP, C1)

        # Compute DHDZ and DZ
        DHDZ[JJ] = -(H4 * DPA * AP * DERIV[JJ] + FC * (H[JJ] - HF)) / G
        DZ = -H4 / (ENMX2 * DHDZ[JJ])

        H_new = H[JJ-1] + DHDZ[JJ] * DZ

        # Overshoot handling
        if H_new >= HV:
            DZ = (HV - H[JJ-1]) / DHDZ[JJ]
            H[JJ] = HV
            Z[JJ] = Z[JJ - 1] + DZ
            WFV = (H[JJ] - HL) / (HV - HL)
            print_lqvp_data(Z[JJ], TEMP, H[JJ], WFV)
            break

        # Normal update
        H[JJ] = H_new
        Z[JJ] = Z[JJ - 1] + DZ
        WFV = (H[JJ] - HL) / (HV - HL)
        print_lqvp_data(Z[JJ], TEMP, H[JJ], WFV)

        JJ += 1

    # Trim the output arrays of the trailing zeros
    Z_trim          = Z[:np.nonzero(Z)[0][-1] + 1]
    T_trim          = T[:np.nonzero(T)[0][-1] + 1]
    P_trim          = P[:np.nonzero(P)[0][-1] + 1]
    H_trim          = H[:np.nonzero(H)[0][-1] + 1]
    DHDZ_trim       = DHDZ[:np.nonzero(DHDZ)[0][-1] + 1]
    DERIV_trim      = DERIV[:np.nonzero(DHDZ)[0][-1] + 1]

    return Z_trim, T_trim, P_trim, H_trim, DHDZ_trim, DERIV_trim


















