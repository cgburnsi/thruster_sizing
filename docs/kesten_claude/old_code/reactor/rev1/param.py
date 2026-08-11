
def PARAM(T, ZA, LOP, CC, HR, LVOP, G, GMMA, K, DPA, BETA, FC, AGM, BGM):
    if (ZA-Z0) < 0.0:
        G = G0 + FC*ZA
    else:
        G = GATZ0
        FC = 0.0
        
    if LVOP == 1:
        GMMA = BGM/T
        # Calculate K, DPA for NH3
        K = AKPHA2*np.exp(-GMMA)/C1**1.6
        DPV=DIF3*(T/492.0)**1.832 * (14.7/PRES)*(1.0-np.exp(-0.0672*(PRES*492.0)/(14.7*T)))
        DPA=DPV
        BETA=(-CC*HR*DPA)/(KP*T)
    else:
        GMMA = BGM/T
        # Calculate K, DPA for N2H4
        K = AKPHA1*np.exp(-GMMA)
        DPL=DIF4*(T/492.0)**1.832 * (14.7/PRES)*(1.0-np.exp(-0.0672*(PRES*492.0)/(14.7*T)))
        DPA=DPL
        BETA=(-CC*HR*DPA)/(KP*T)

