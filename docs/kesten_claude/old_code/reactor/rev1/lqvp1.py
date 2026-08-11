from unbar import UNBAR
from param01 import PARAM
from print_data import print_lqvp_data

def LQVP(H, ZLV, Q, JJ, Q1, TEMP,
         DERIV, DHDZ, Z,
         FC, HF, HL, HV, ENMX2, G0, TVAP, Z0, DZ_init,
         AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1,
         max_steps=250):

    # First entry: use ZLV and input values
    DERIV[JJ] = Q
    DHDZ[JJ] = Q1
    Z[JJ] = ZLV
    JJ += 1

    while JJ < max_steps:
        TEMP = TVAP
        DERIV[JJ] = DERIV[JJ - 1]

        # Use previous Z for lookups
        Z_prev = Z[JJ - 1]
        H4 = UNBAR('TBLH4', TEMP)
        AP = UNBAR('ZTBLAP', Z_prev)
        A  = UNBAR('ZTBLA', Z_prev)

        GATZ0 = G0 + FC * Z0
        G, GMMA, K, BETA, DPA = PARAM(
            TEMP, Z_prev, 1, 0.0, 0.0, 0,
            Z0, G0, FC, GATZ0,
            AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4,
            PRES, KP, C1
        )

        # Compute DHDZ and DZ
        DHDZ[JJ] = -(H4 * DPA * AP * DERIV[JJ] + FC * (H - HF)) / G
        DZ = -H4 / (ENMX2 * DHDZ[JJ])

        H_new = H + DHDZ[JJ] * DZ

        # Overshoot handling
        if H_new >= HV:
            DZ = (HV - H) / DHDZ[JJ]
            H = HV
            Z[JJ] = Z[JJ - 1] + DZ
            WFV = (H - HL) / (HV - HL)
            print_lqvp_data(Z[JJ], TEMP, H, WFV)
            break

        # Normal update
        H = H_new
        Z[JJ] = Z[JJ - 1] + DZ
        WFV = (H - HL) / (HV - HL)
        print_lqvp_data(Z[JJ], TEMP, H, WFV)

        JJ += 1

    return H, JJ, TEMP
