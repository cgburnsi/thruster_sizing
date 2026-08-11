from print_data import print_inputs, print_liquid_data

import numpy as np
from unbar import UNBAR
from param2 import PARAM


# Constants
FRAC        = 0.99      # Tolerance Factor in trapezoidal calculation
TOL         = 0.01      # Sets a convergence criteria margin around the target value                  


def LIQUID(state, config):
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
    ENMX1   = config.ENMX1
    TF      = config.TF
    CFL     = config.CFL
    
    
    # Retrieve State Parameters
    Z       = state.Z
    H       = state.H
    T       = state.TEMP
    P       = state.PRES
    DHDZ    = state.DHDZ
    DERIV   = state.DERIV
    
    DZ          = 0.0
    Z[0]        = config.Z0
    H[0]        = HF
    T[0]        = TF
    P[0]        = config.PRES
    DHDZ[0]     = 0.0
    GATZ0       = G0 + FC * Z0
    
    TVAP        = UNBAR('TVAP', P[0])
    HL          = (TVAP - TF) * CFL
    VP          = UNBAR('TBLVP', T[0])
    CN2H4       = (VP * config.WM4) / (config.R * T[0])
    H4          = UNBAR('TBLH4', T[0])
    AP          = UNBAR('ZTBLAP', Z[0])
    A           = UNBAR('ZTBLA', Z[0])
    
    II = 1
    while II < config.LIQ_MAX:
    
        Z[II] = Z[II - 1] + DZ
        T[II] = TF + (H[II-1] - HF) / CFL
    
        P[II] = P[II - 1]
    
    
        VP    = UNBAR('TBLVP', T[II])
        CN2H4 = (VP * WM4) / (R * T[II])
        H4    = UNBAR('TBLH4', T[II])
        AP    = UNBAR('ZTBLAP', Z[II])
        A     = UNBAR('ZTBLA', Z[II])
    
        G, GMMA, K, BETA, DPA = PARAM(T[II], Z[II], 1, CN2H4, H4, 0, Z0, G0, FC, GATZ0, AGM, BGM, ALPHA1, ALPHA2, DIF3, DIF4, P[II], KP, C1)
        DERIV[II], MI = SLOPE(CN2H4, GMMA, K, BETA, DPA, A, EN1, HL, DIF3, T[II], P[II], G, AP, WM4, R)
    
        #print_liquid_data(Z[II], T[II], H[II-1], DHDZ[II])
    
    
        if MI >= 20:
            print('MI >= 20')
            DERIV[II] = DERIV[II - 1]
    
        DHDZ[II] = -(H4 * DPA * AP * DERIV[II] + FC * (H[II] - HF)) / G
        DZ = -H4 / (ENMX1 * DHDZ[II])
    
        # Predict next enthalpy
        H_next = H[II - 1] + DHDZ[II] * DZ
        
        # Overshoot correction
        if H_next >= HL:
            DZ += (HL - H_next) / DHDZ[II]
            H_next = HL
        H[II] = H_next
    
        print_liquid_data(Z[II], T[II], H[II-1], DHDZ[II])
        
        if H[II-1] >= HL:
            break
    
        II += 1
        
    # Trim the output arrays of the trailing zeros
    Z_trim          = Z[:np.nonzero(Z)[0][-1] + 1]
    T_trim          = T[:np.nonzero(T)[0][-1] + 1]
    P_trim          = P[:np.nonzero(P)[0][-1] + 1]
    H_trim          = H[:np.nonzero(H)[0][-1] + 1]
    DHDZ_trim       = DHDZ[:np.nonzero(DHDZ)[0][-1] + 1]
    DERIV_trim      = DERIV[:np.nonzero(DHDZ)[0][-1] + 1]

    return Z_trim, T_trim, P_trim, H_trim, DHDZ_trim, DERIV_trim








def get_integration_indices(MM, ADIV:int = 100, BDIV:int = 100):
    ''' Provides the appropriate indexing parameters based on the 
        current MM iteration value. MMAX is the total number of itegration
        points.  MINT is the index where the second region [B, A] begins.
        MAX is the last valid index. '''
    if MM == 0:
        MINT = 0
        MMAX = BDIV + 1
    else:
        MINT = ADIV + 1
        MMAX = ADIV + BDIV + 1

    MAX = MMAX - 1
    
    return MINT, MMAX, MAX

def has_converged(DERIF, prev_deriv, IK, KJ, tol=0.05):
    if KJ == 1 and prev_deriv is not None:
        return abs(DERIF - prev_deriv) <= tol * DERIF and IK == 1
    return False

def check_convergence(I, H, HL, CG, C4, DIFF, TEMP, PRES, G, AP, DPA, 
                      WMAV, R, KK, DERIF, PRINT, XOA, CPA, A, B):
    """
    Check convergence of surface reaction rate in the SLOPE function.

    Parameters:
    - I: index into CG, CPA, DPA, etc.
    - H: current enthalpy
    - HL: latent heat boundary enthalpy
    - CG: concentration array
    - C4: some threshold or target concentration
    - DIFF: base diffusion coefficient
    - TEMP: temperature in degR
    - PRES: pressure in psia
    - G: mass flux
    - AP: particle surface area
    - DPA: pressure difference array
    - WMAV: average molecular weight
    - R: gas constant
    - VISVST: viscosity table
    - KK: table size or additional context for viscosity
    - DERIF: surface gradient computed from internal profile
    - PRINT: flag for debug printing
    - XOA, CPA: surface mole fraction and concentration profiles
    - A, B: arrays or scalars related to geometry

    Returns:
    - RATE: updated rate
    - MFLAG: convergence flag (1 = converged, 0 = not converged)
    """
    RATE = DERIF

    if C4 != CG:

        return RATE, None  # Skip check

    if H <= HL:

        return RATE, None  # Still in liquid region, skip check

    # Adjusted diffusion coefficient (empirical law)
    DIFN = DIFF * ((TEMP / 492.0) ** 1.823 * 14.7 / PRES)

    # Get viscosity at current TEMP
    VIS = UNBAR('VISVST', TEMP)

    # Compute density from ideal gas law
    RHO = PRES * WMAV / (R * TEMP)

    # Chilton-Colburn analogy for mass transfer coefficient
    AKC = (0.61 * G) / RHO * ( (VIS / (RHO * DIFN)) ** -0.667 ) * ( (G / (AP * VIS)) ** -0.41 )

    # Compute external rate estimate
    RAT = AKC * CG / DPA

    if RATE > RAT:
        MFLAG = 1
        RATE = RAT
    else:
        MFLAG = 0


    # Optional diagnostic output
    if B > 0.998 * A:
        if PRINT == 1:
            print(f"ITERATION = {I}")
            print("        X/A       CPA        X/A       CPA       X/A       CPA       X/A       CPA")
            for i in range(0, len(XOA), 4):
                segment = "  ".join([f"{XOA[j]:12.7e} {CPA[j]:12.7e}" for j in range(i, min(i+4, len(XOA)))])
                print("  ", segment)
            print(f"\n        THE SLOPE CONVERGES TO {RATE:12.7e}")

    return RATE, MFLAG







def compute_reaction_profile(CG, GMMA, K, BETA, EN12, A, B, DR, BDR, MINT, MMAX):
    ''' Computes and returns the weighted local reaction rate profile S[i] 
        across the catalyst particle radius. Fully functional — no mutation. '''
    
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
    ''' Performs a trapezoidal integration of the reaction rate profile over the radial dimension 
        of a porous catalyst particle.  It splits the integration domain into two sections [0, B]
        and [B, A].  This split allows the trapezoidal rule to apply correct weights per region.  
        SUM is the summation of S in region [0, B], and SUMA is in the region [B, A]. 
        The loop indexing starts at 1.  The trapezoidal rule doesn't use the 0 index in the
        calculation.  
        
        D2 is used in other functions to determine if the penetration depth B is good enough.
        TOT is is used to softe the convergence criteria to prevent oscillations and slow convergence. '''
        
    SUM     = 0.0                   # Accumulated integral over the region [0, B]
    SUMA    = 0.0                   # Accumulated integral over the region [B, A]
    S[0]    = 0.0                   # Weighted Local Reaction Rate is not needed in the calculation
    
    SUM  = sum(S[1:MINT+1])
    SUMA = sum(S[MINT+1:MAX])
    
    # Compute SGMA
    sgma_first  = (S[1] + 2.0 * SUM) * (DR / (2.0 * DPA))       
    sgma_second = (S[MMAX-1] + 2.0 * SUMA) * (BDR / (2.0 * DPA))  
    SGMA        = sgma_first + sgma_second
    D2          = SGMA - FRAC * CG
    TOT         = D2 + TOL * CG
    
    return D2, TOT



def SLOPE(CG, GMMA, K, BETA, DPA, A, EN12, HL, DIFF3, TEMP, PRES, G0, AP, WM4, R):

    ADIV, BDIV = 100, 100                               # Number of points in each region [0, B] and [B, A]

    # CG is the concentration of reactant gas at the particle surface    


    # Local arrays sized to the maximum expected number of points.  In
    # the original Fortran code these were dimensioned to 210 or 250
    # depending on usage.  Here we allocate a generous upper bound (300)
    # and will slice as needed.
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

    # Flags and counters
    IFLAG = 0
    JFLAG = 0

    # Starting guess for B and loop counters
    B = 0.0
    I = 1
    
    STORE_VAL = 1.0
    KJ = 0
    HOLD = 0.0
    NI = 0
    MI = 0
    IK = 0
    MM = 0
    IM = 0
    F_factor = 0.5

    # Variables that are assigned conditionally in the Fortran code
    # but referenced in general expressions must be initialised.
    # XOLD stores the previous B value when adjusting B.
    XOLD = 0.0
    # D1 stores the previous D2 value when adjusting B.
    D1 = 0.0
    

    # Main iteration loop.  This corresponds to labels 20, 110, 115,
    # etc. in the original Fortran routine.  The loop continues until
    # the derivative converges or a maximum number of iterations is
    # reached.
    RATE = 0.0
    # Save the previous derivative to check convergence
    prev_deriv = None
    # Set an upper bound on the total number of iterations to avoid
    # infinite loops in the case of non‑convergence.
    max_total_iters = 200

    total_iters = 0
        
    while True:
        total_iters += 1
        if total_iters > max_total_iters:
            # If convergence has not been achieved, return the past computed RATE and break out of the loop
            return RATE, MI
                # Increment iteration counter
        MM += 1
        
        # Step 3: Domain Discritization
        DR              = B / ADIV                                      #  Discritization of domain from 0 to B (particle center to B)
        BDR             = (A - B) / BDIV                                #  Discritization of domain from B to A (slope change to particle surface)

        
        MINT, MMAX, MAX = get_integration_indices(MM, ADIV, BDIV)
        S               = compute_reaction_profile(CG, GMMA, K, BETA, EN12, A, B, DR, BDR, MINT, MMAX)
        D2, TOT         = integrate_reaction_profile(CG, S, DR, BDR, MINT, MAX, MMAX, DPA)              # Trapezoidal Integration
        
        # Handle early iterations and D2-based checks
        if MM <= 1:
            #print('MM <=1')
            if IFLAG == 1:
                #print('MM <= 1, IFLAG==1')
                pass  # continue below (label 115 in original)
            elif D2 < 0:
                #print('MM <= 1, D2 < 0')
                pass  # original label 150, handled later
            # otherwise fall through to main flow

        # Decide whether to proceed to 115 or branch to simplified model (230)
        if D2 <= 0.0 and TOT > 0.0:
            #print('D2 <= 0.0 and TOT > 0.0')
            branch_to_230 = True
        else:
            #print('else 1')
            # Continue iterating
            if MM <= 0:
                #print('MM <= 0')
                continue  # just go to next iteration

            if MM == 0:
                #print('MM == 0')
                continue  # Fortran label 120

            if MM - 1 == 0:
                #print('MM - 1 == 0')
                # Save old values and nudge B closer to A
                XOLD = B
                D1 = D2
                B = 0.999999 * A
                continue  # restart iteration
            else:
                #print('else 2')
                # Improve estimate for B using linear extrapolation
                TJMP = B
                MI += 1
                if (D2 - D1) != 0:
                    B += (D2 / (D2 - D1)) * (XOLD - B)
                XOLD = TJMP
                D1 = D2

                if MI > 20:
                    #print('MI > 20')
                    if JFLAG == 1:
                        #print('MI > 20, JFLAG ==1')
                        branch_to_230 = True
                    else:
                        #print('else 3')
                        B = 0.9 * A
                        JFLAG = 1
                        # Reset state
                        MM = 0
                        KJ = 0
                        IFLAG = 0
                        continue  # restart iteration
                else:
                    #print('else 4')
                    continue  # restart iteration

        
      
        
      
        
      
        
        
        # The following code is what I kinda say is the second part of the SLOPE function.
        
        
        if branch_to_230:
            #print('branch_to_230')
            # Equivalent to label 230: report that a satisfactory starting
            # curve has been found and compute derivative with simplified
            # approach
            # In Fortran several I/O operations occur here; we omit them.
            if B > 0.998 * A:
                #print('B > 0.998 * A')
                DERIF = 2.6 * CG / (A - B) if (A - B) != 0 else 0.0
                RATE = DERIF

                # Call convergence check against external rate estimate
                RATE, MFLAG = check_convergence(I=I, H=HOLD, HL=HL,
                    CG=CG,
                    C4=CG,
                    DIFF=DIFF3,
                    TEMP=TEMP,
                    PRES=PRES,
                    G=G0,
                    AP=AP,
                    DPA=DPA,
                    WMAV=WM4,
                    R=R,
                    KK=70,
                    DERIF=DERIF,
                    PRINT=1,
                    XOA=XOA,
                    CPA=CPA,
                    A=A,
                    B=B
                )


                return RATE, MI
            # Continue with range of X from zero to A (label 237)
            X_val = 0.0
            NI = 1
            # Compute CPB and other arrays for the full range 0..A
            for NN in range(1, MMAX + 1):
                print('NN in range(1, MMAX+1')
                if NN > MINT:
                    KJK = NN - MINT
                    X_val = B + float(KJK) * BDR
                else:
                    X_val = float(NN - 1) * DR
                # When IM >= 1 use CPA from previous iteration, otherwise
                # compute CPB directly
                if IM >= 1:
                    print('IM >= 1')
                    pass  # skip updating CPB
                else:
                    print('IM >= 1 ----- ELSE')
                    denom = A - B
                    if denom != 0:
                        CPB[NN] = CG * X_val / denom - (B / denom) * CG
                    else:
                        CPB[NN] = 0.0
                    if CPB[NN] < 0.0:
                        CPB[NN] = 0.0
                # Evaluate TERM1, TERM2, FST, SEC for CPB (or CPA if IM>=1)
                cp_val = CPA[NN] if IM >= 1 else CPB[NN]
                TERM1[NN] = K * (CG ** (1.0 - EN12)) * (cp_val ** EN12)
                beta_val = BETA
                gmma_val = GMMA
                if CG != 0.0:
                    print('CG != 0.0')
                    exponent = gmma_val * beta_val * (1.0 - cp_val / CG) / (1.0 + beta_val * (1.0 - cp_val / CG))
                else:
                    print('CG != 0.0 ----- ELSE')
                    exponent = 0.0
                TERM2[NN] = np.exp(exponent)
                XXX[NN] = X_val
                XOA[NN] = XXX[NN] / A if A != 0 else 0.0
                RCP = TERM1[NN] * TERM2[NN]
                FST[NN] = X_val * X_val * RCP
                if X_val != 0.0:
                    SEC[NN] = (1.0 / X_val - 1.0 / A) * X_val * X_val * RCP
                else:
                    SEC[NN] = 0.0
            # After computing CPB (and possibly CPA) arrays we fall
            # through to the trapezoidal integration at label 172
        # If we reach here without branching we are in the main iteration
        # path that computes integrals using CPA from the previous
        # iteration when KJ==1 or IM>0.
        # Compute TERM1, TERM2, FST, SEC for all points using CPA
        if KJ == 1:
            print('KJ == 1')
            # At label 150 in Fortran the initial guess through the
            # origin is satisfactory; compute the integrals for
            # CPA=0 at X=0 and then iterate.
            # Here we reuse the values already computed above when
            # branch_to_230 was False and KJ==1.
            pass
        # Compute trapezoidal integrals C and D
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
        # Compute E and CPA arrays
        E[0] = D[MMAX-1]
        CPA[0] = 0.0
        STORE_VAL = CPA[0]
        if not (KJ == 1 or IM == 0 or NI == 1):
            pass  # skip CPA update
        SAM = 0.0
        for LL in range(1, MMAX):
            print('LL in range(1, MMAX)')
            if LL > MINT:
                E[LL] = E[LL - 1] - (SEC[LL] + SEC[LL - 1]) * (BDR / (2.0 * DPA)) if DPA != 0 else 0.0
                SAM += BDR
            else:
                E[LL] = E[LL - 1] - (SEC[LL] + SEC[LL - 1]) * (DR / (2.0 * DPA)) if DPA != 0 else 0.0
                SAM = float(LL - 1) * DR
            # Compute new CPA
            if SAM != 0.0:
                CPA[LL] = CG - ((1.0 / SAM) - (1.0 / A)) * C[LL] - E[LL]
            else:
                CPA[LL] = 0.0
            if LL < MINT:
                CPA[LL] = 0.0
            if CPA[LL] < 0.0:
                CPA[LL] = 0.0
        # After computing CPA decide next iteration path
        if KJ == 1:
            # go to 280
            pass
        # Determine whether to blend CPA with CPB (label 250)
        if IM == 0:
            # skip blending, go to label 250 equivalent
            pass
        else:
            IK = 1
        # If IK == 1 jump to label 280 for derivative calculation
        # If not, prepare to recalculate TERM arrays using CPA (label 192)
        if IK != 1:
            # recompute TERM arrays using updated CPA
            for LI in range(MMAX):
                if LI > MINT:
                    KLK = LI - MINT
                    X_val = B + float(KLK) * BDR
                else:
                    X_val = float(LI - 1) * DR
                TERM1[LI] = K * (CG ** (1.0 - EN12)) * (CPA[LI] ** EN12)
                beta_val = BETA
                gmma_val = GMMA
                if CG != 0.0:
                    exponent = gmma_val * beta_val * (1.0 - CPA[LI] / CG) / (1.0 + beta_val * (1.0 - CPA[LI] / CG))
                else:
                    exponent = 0.0
                TERM2[LI] = np.exp(exponent)
                XXX[LI] = X_val
                XOA[LI] = XXX[LI] / A if A != 0 else 0.0
                RCP = TERM1[LI] * TERM2[LI]
                FST[LI] = X_val * X_val * RCP
                if X_val != 0.0:
                    SEC[LI] = (1.0 / X_val - 1.0 / A) * X_val * X_val * RCP
                else:
                    SEC[LI] = 0.0
            # Return to trapezoidal integration loop
            continue
        # Blend CPA with CPB when appropriate (label 250)
        for IL_idx in range(MMAX):
            CPA[IL_idx] = F_factor * CPA[IL_idx] + (1.0 - F_factor) * CPB[IL_idx]
        IK = 0
        # Compute derivative at point near A(I)
        if MMAX - 1 >= 1:
            DERIF = (CG - CPA[MAX]) / BDR if BDR != 0 else 0.0
        else:
            DERIF = 0.0
        IM += 1
        if IM > 99 and IFLAG == 1:
            # simplified version does not converge; reset B and start over
            B = 0.000001 * A
            DR = B / ADIV if ADIV != 0 else 0.0
            BDR = (A - B) / BDIV if BDIV != 0 else 0.0
            MMAX = ADIV + BDIV + 1
            MINT = ADIV + 1
            MAX = MMAX - 1
            IM = 0
            KJ = 0
            # Jump to label 237
            continue
        elif IM > 99:
            # if not converged in 99 iterations, continue with convergence check
            pass
        # Check convergence of derivative
        if KJ == 1:
            IK = 1
        if prev_deriv is not None:
            if abs(DERIF - prev_deriv) <= 0.05 * DERIF:
                if IK == 1:
                    RATE = DERIF
                    break
        prev_deriv = DERIF
        if KJ == 1:
            # Jump to label 192
            continue
        if IK == 1:
            # Jump to label 250
            continue
        # Default: jump to label 192
        continue
    # End of while loop
    # Additional checks based on Fortran logic for C4 and HL have been
    # omitted.  The final RATE value is returned.
    return RATE, MI
        
        
        
        
        
        
        
        
        
        

    
    
    
    