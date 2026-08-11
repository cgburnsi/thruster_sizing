import numpy as np
from config import ReactorConfig
from properties import PropertyTables

# Constants
FRAC = 0.99
TOL = 0.01

def get_integration_indices(MM, ADIV:int = 100, BDIV:int = 100):
    if MM == 0:
        MINT = 0
        MMAX = BDIV + 1
    else:
        MINT = ADIV + 1
        MMAX = ADIV + BDIV + 1
    MAX = MMAX - 1
    return MINT, MMAX, MAX

def check_convergence(I, H, HL, CG, C4, DIFF, TEMP, PRES, G, AP, DPA, 
                      WMAV, R, KK, DERIF, PRINT, XOA, CPA, A, B, props: PropertyTables):
    """
    Check convergence of surface reaction rate in the SLOPE function.
    """
    RATE = DERIF

    if C4 != CG:
        return RATE, None

    if H <= HL:
        return RATE, None

    # Adjusted diffusion coefficient
    DIFN = DIFF * ((TEMP / 492.0) ** 1.823 * 14.7 / PRES)

    # Viscosity
    VIS = props.VISVST(TEMP)

    # Density
    RHO = PRES * WMAV / (R * TEMP)

    # Safety checks
    if RHO == 0: RHO = 1e-30
    if VIS == 0: VIS = 1e-30
    if AP == 0: AP = 1e-30
    if DIFN == 0: DIFN = 1e-30

    # Mass transfer coefficient (Chilton-Colburn)
    AKC = (0.61 * G) / RHO * ( (VIS / (RHO * DIFN)) ** -0.667 ) * ( (G / (AP * VIS)) ** -0.41 )

    # External rate estimate
    if DPA == 0:
        RAT = 0.0
    else:
        RAT = AKC * CG / DPA

    if RATE > RAT:
        MFLAG = 1
        RATE = RAT
    else:
        MFLAG = 0

    return RATE, MFLAG


def compute_reaction_profile(CG, GMMA, K, BETA, EN12, A, B, DR, BDR, MINT, MMAX):
    S = [0.0] * MMAX

    for IJ in range(MMAX):
        if IJ > MINT:
            KMJ = IJ - MINT
            X_val = B + KMJ * BDR
        else:
            X_val = (IJ - 1) * DR if IJ > 0 else 0.0

        denom = A - B
        if denom == 0.0:
            CP = 0.0
        else:
            CP = CG * X_val / denom - (B / denom) * CG
        CP = max(CP, 0.0)

        FACT1 = K * (CG ** (1.0 - EN12)) * (CP ** EN12)
        if CG != 0.0:
            exponent = GMMA * BETA * (1.0 - CP / CG) / (1.0 + BETA * (1.0 - CP / CG))
        else:
            exponent = 0.0
        FACT2 = np.exp(exponent)

        if X_val != 0.0:
            RCP = FACT1 * FACT2
            S[IJ] = (1.0 / X_val - 1.0 / A) * X_val**2 * RCP
        else:
            S[IJ] = 0.0

    return S

def integrate_reaction_profile(CG, S, DR, BDR, MINT, MAX, MMAX, DPA):
    if DPA == 0: DPA = 1e-30
    
    SUM  = sum(S[1:MINT+1])
    SUMA = sum(S[MINT+1:MAX])
    
    sgma_first  = (S[1] + 2.0 * SUM) * (DR / (2.0 * DPA))       
    sgma_second = (S[MMAX-1] + 2.0 * SUMA) * (BDR / (2.0 * DPA))  
    SGMA        = sgma_first + sgma_second
    D2          = SGMA - FRAC * CG
    TOT         = D2 + TOL * CG
    
    return D2, TOT


def SLOPE(CG, GMMA, K, BETA, DPA, A, TEMP, G, AP, cfg: ReactorConfig, props: PropertyTables):
    """
    Solves the internal diffusion-reaction problem for the liquid phase.
    Refactored to use ReactorConfig.
    """
    
    # Unpack config constants needed for calculation
    EN12 = cfg.EN1
    HL = cfg.HL
    DIFF3 = cfg.DIF3
    PRES = cfg.PRES
    WM4 = cfg.WM4
    R = cfg.R
    
    ADIV, BDIV = 100, 100
    
    MAX_SIZE = 300
    FST = [0.0] * MAX_SIZE
    SEC = [0.0] * MAX_SIZE
    D = [0.0] * MAX_SIZE
    CPA = [0.0] * MAX_SIZE
    CPB = [0.0] * MAX_SIZE
    E = [0.0] * MAX_SIZE
    TERM1 = [0.0] * MAX_SIZE
    TERM2 = [0.0] * MAX_SIZE
    C = [0.0] * MAX_SIZE
    XXX = [0.0] * MAX_SIZE
    XOA = [0.0] * MAX_SIZE

    IFLAG = 0
    JFLAG = 0
    B = 0.0
    I = 1
    KJ = 0
    HOLD = 0.0
    NI = 0
    MI = 0
    IK = 0
    MM = 0
    IM = 0
    F_factor = 0.5
    XOLD = 0.0
    D1 = 0.0
    
    RATE = 0.0
    prev_deriv = None
    max_total_iters = 200
    total_iters = 0
        
    while True:
        total_iters += 1
        if total_iters > max_total_iters:
            return RATE, MI

        MM += 1
        
        DR = B / ADIV                                      
        BDR = (A - B) / BDIV                                
        
        MINT, MMAX, MAX = get_integration_indices(MM, ADIV, BDIV)
        S = compute_reaction_profile(CG, GMMA, K, BETA, EN12, A, B, DR, BDR, MINT, MMAX)
        D2, TOT = integrate_reaction_profile(CG, S, DR, BDR, MINT, MAX, MMAX, DPA)
        
        # Logic for B adjustment
        branch_to_230 = False
        if D2 <= 0.0 and TOT > 0.0:
            branch_to_230 = True
        else:
            if MM <= 0: continue 
            if MM == 0: continue 

            if MM - 1 == 0:
                XOLD = B
                D1 = D2
                B = 0.999999 * A
                continue 
            else:
                TJMP = B
                MI += 1
                if (D2 - D1) != 0:
                    B += (D2 / (D2 - D1)) * (XOLD - B)
                XOLD = TJMP
                D1 = D2

                if MI > 20:
                    if JFLAG == 1:
                        branch_to_230 = True
                    else:
                        B = 0.9 * A
                        JFLAG = 1
                        MM = 0
                        KJ = 0
                        IFLAG = 0
                        continue  
                else:
                    continue  

        if branch_to_230:
            if B > 0.998 * A:
                DERIF = 2.6 * CG / (A - B) if (A - B) != 0 else 0.0
                RATE = DERIF

                # Convergence check
                RATE, MFLAG = check_convergence(
                    I=I, H=HOLD, HL=HL,
                    CG=CG, C4=CG, DIFF=DIFF3, TEMP=TEMP, PRES=PRES,
                    G=G, AP=AP, DPA=DPA, WMAV=WM4, R=R, KK=70,
                    DERIF=DERIF, PRINT=1, XOA=XOA, CPA=CPA, A=A, B=B,
                    props=props
                )
                return RATE, MI

            # Calculate Profile
            X_val = 0.0
            NI = 1
            for NN in range(1, MMAX + 1):
                if NN > MINT:
                    KJK = NN - MINT
                    X_val = B + float(KJK) * BDR
                else:
                    X_val = float(NN - 1) * DR
                
                denom = A - B
                if denom != 0:
                    CPB[NN] = CG * X_val / denom - (B / denom) * CG
                else:
                    CPB[NN] = 0.0
                if CPB[NN] < 0.0: CPB[NN] = 0.0
                
                cp_val = CPA[NN] if IM >= 1 else CPB[NN]
                
                # Reaction Terms
                TERM1[NN] = K * (CG ** (1.0 - EN12)) * (cp_val ** EN12)
                
                if CG != 0.0:
                    exponent = GMMA * BETA * (1.0 - cp_val / CG) / (1.0 + BETA * (1.0 - cp_val / CG))
                else:
                    exponent = 0.0
                TERM2[NN] = np.exp(exponent)
                
                XXX[NN] = X_val
                RCP = TERM1[NN] * TERM2[NN]
                FST[NN] = X_val * X_val * RCP
                if X_val != 0.0:
                    SEC[NN] = (1.0 / X_val - 1.0 / A) * X_val * X_val * RCP
                else:
                    SEC[NN] = 0.0
            
        C[0] = 0.0
        for JJ in range(1, MMAX):
            if JJ > MINT:
                C[JJ] = C[JJ - 1] + (FST[JJ] + FST[JJ - 1]) * (BDR / (2.0 * DPA)) if DPA != 0 else 0.0
            else:
                C[JJ] = C[JJ - 1] + (FST[JJ] + FST[JJ - 1]) * (DR / (2.0 * DPA)) if DPA != 0 else 0.0
        
        SEC[0] = 0.0
        D[0] = 0.0
        for KK in range(1, MMAX):
            if KK > MINT:
                D[KK] = D[KK - 1] + (SEC[KK - 1] + SEC[KK]) * (BDR / (2.0 * DPA)) if DPA != 0 else 0.0
            else:
                D[KK] = D[KK - 1] + (SEC[KK - 1] + SEC[KK]) * (DR / (2.0 * DPA)) if DPA != 0 else 0.0
        
        E[0] = D[MMAX-1]
        CPA[0] = 0.0
        SAM = 0.0
        for LL in range(1, MMAX):
            if LL > MINT:
                E[LL] = E[LL - 1] - (SEC[LL] + SEC[LL - 1]) * (BDR / (2.0 * DPA)) if DPA != 0 else 0.0
                SAM += BDR
            else:
                E[LL] = E[LL - 1] - (SEC[LL] + SEC[LL - 1]) * (DR / (2.0 * DPA)) if DPA != 0 else 0.0
                SAM = float(LL - 1) * DR
            
            if SAM != 0.0:
                CPA[LL] = CG - ((1.0 / SAM) - (1.0 / A)) * C[LL] - E[LL]
            else:
                CPA[LL] = 0.0
            if LL < MINT:
                CPA[LL] = 0.0
            if CPA[LL] < 0.0:
                CPA[LL] = 0.0
        
        if IM != 0: IK = 1
        
        if IK != 1:
            # Recompute FST/SEC with new CPA
            for LI in range(MMAX):
                if LI > MINT:
                    KLK = LI - MINT
                    X_val = B + float(KLK) * BDR
                else:
                    X_val = float(LI - 1) * DR
                    
                TERM1[LI] = K * (CG ** (1.0 - EN12)) * (CPA[LI] ** EN12)
                if CG != 0.0:
                    exponent = GMMA * BETA * (1.0 - CPA[LI] / CG) / (1.0 + BETA * (1.0 - CPA[LI] / CG))
                else:
                    exponent = 0.0
                TERM2[LI] = np.exp(exponent)
                RCP = TERM1[LI] * TERM2[LI]
                FST[LI] = X_val * X_val * RCP
                if X_val != 0.0:
                    SEC[LI] = (1.0 / X_val - 1.0 / A) * X_val * X_val * RCP
                else:
                    SEC[LI] = 0.0
            continue
        
        for IL_idx in range(MMAX):
            CPA[IL_idx] = F_factor * CPA[IL_idx] + (1.0 - F_factor) * CPB[IL_idx]
        IK = 0
        
        if MMAX - 1 >= 1:
            DERIF = (CG - CPA[MAX]) / BDR if BDR != 0 else 0.0
        else:
            DERIF = 0.0
        IM += 1
        
        if IM > 99 and IFLAG == 1:
            B = 0.000001 * A
            DR = B / ADIV if ADIV != 0 else 0.0
            BDR = (A - B) / BDIV if BDIV != 0 else 0.0
            MMAX = ADIV + BDIV + 1
            MINT = ADIV + 1
            MAX = MMAX - 1
            IM = 0
            KJ = 0
            continue
        
        if KJ == 1: IK = 1
        if prev_deriv is not None:
            if abs(DERIF - prev_deriv) <= 0.05 * DERIF:
                if IK == 1:
                    RATE = DERIF
                    break
        prev_deriv = DERIF
        if KJ == 1: continue
        if IK == 1: continue
        continue
    
    return RATE, MI