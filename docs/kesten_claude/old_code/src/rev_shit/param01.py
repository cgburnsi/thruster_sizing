import numpy as np

def PARAM(T, ZA, LOP, CC, HR, LVOP, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, PRES, KP, C1):

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

    # The comment questions if it should be CC+HR, not CC*HR. We'll keep CC+HR for accuracy to code.
    BETA = -(CC * HR * DPA) / (KP * T)

    return G, GMMA, K, BETA, DPA