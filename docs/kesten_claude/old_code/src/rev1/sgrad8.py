import numpy as np
from unbar import UNBAR


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

    print("\n🔍 Entering TRAPP integration loop...")
    print(f"{'i':>3} {'X':>10} {'CPX':>10} {'RHET':>10} {'FOXI1':>10}")

    # First point
    CPX1 = CPXF(U, XOA, CPS)
    RHET1 = RHETF(CPS, CPX1, GAMMA, BETA, K0, 1)
    TRM1 = FOXI1(U, RHET1) / 2.0
    print(f"{0:>3} {U:10.5e} {CPX1:10.5e} {RHET1:10.5e} {FOXI1(U, RHET1):10.5e}")

    # Internal points
    xs = np.linspace(U + H, V - H, N - 1)
    for i, x in enumerate(xs, start=1):
        CPX = CPXF(x, XOA, CPS)
        RHET = RHETF(CPS, CPX, GAMMA, BETA, K0, 1)
        FX = FOXI1(x, RHET)
        print(f"{i:>3} {x:10.5e} {CPX:10.5e} {RHET:10.5e} {FX:10.5e}")
        SUM += FX

    # Last point
    CPX2 = CPXF(V, XOA, CPS)
    RHET2 = RHETF(CPS, CPX2, GAMMA, BETA, K0, 1)
    TRM2 = FOXI1(V, RHET2) / 2.0
    print(f"{N:>3} {V:10.5e} {CPX2:10.5e} {RHET2:10.5e} {FOXI1(V, RHET2):10.5e}")

    RIESUM = H * (TRM1 + SUM + TRM2)

    print(f"\n📏 TRAPP complete. H = {H:.5e}, RIESUM = {RIESUM:.5e}\n")
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

def is_converged(TMTPO, TMTPN, CMCPO, CMCPN):
    dT_rel = abs(TMTPO - TMTPN) / max(TMTPN, 1e-6)
    dC_rel = abs(CMCPO - CMCPN) / max(CMCPN, 1e-6)
    print(f"   ΔT/T = {dT_rel:.3e}")
    print(f"   ΔC/C = {dC_rel:.3e}")
    return dT_rel < 0.05 and dC_rel < 0.05

def temperature_fluctuates(T1, T2, T3):
    return min(T1, T2, T3) < T2 or max(T1, T2, T3) > T2

def SGRAD(TEMP, PRES, GG, C1, C2, C3, C4, DIF3, DIF4, A, AP, BGM, KP, ALPHA2, ALPHA3, EN3):
    MAXITER = 50
    T = TEMP
    P = PRES
    G = GG
    D03, D04 = DIF3, DIF4
    CI1, CI2, CI3, CI4 = C1, C2, C3, C4
    LTFLG = 0
    WAF1 = 0.8
    WAF2 = 1.0 - WAF1
    NPART = 50

    CPS = CI3 / 2.0
    CMCPN = CI3 - CPS
    CMCPO = 0.0
    TMTPO = 0.0
    TMTPN = 0.0
    TPSP = TPSPP = TPS = T
    X0 = 0.0

    while WAF1 <= 0.95:
        print(f"\n🔁 Starting WAF1 = {WAF1:.2f}")
        for LP1 in range(1, MAXITER + 1):
            print(f"\n➡️  Iteration LP1 = {LP1}")
            HC, RHO, MU = CALC_HC(T, G, AP, CI1, CI2, CI3, CI4)
            KC3, KC4 = CALC_KC(T, P, D03, D04, G, RHO, MU, AP)
            DP3 = DP3F(T, D03, P)
            DCPDX = KC3 / DP3 * (CI3 - CPS)
            X0_est = A - CPS / DCPDX
            print(f"   CPS = {CPS:.5e}")
            print(f"   DCPDX = {DCPDX:.5e}")
            print(f"   X0 estimate = {X0_est:.5e}")

            TPS = T - (UNBAR('H4TBL', T) * KC4 * CI4 + UNBAR('H3TBL', T) * DP3 * DCPDX) / HC
            TPS = max(TPS, 1.0)
            H3 = UNBAR('H3TBL', TPS)
            DP3 = DP3F(TPS, D03, P)
            GAMMA = BGM / TPS
            BETA = -CPS * H3 * DP3 / (KP * TPS)
            KO = ALPHA2 * np.exp(-GAMMA) * CI1 ** EN3

            if X0_est < 0:
                print("⚠️  X0 < 0: using fallback CPS and restarting...")
                X0 = 0.0
                CPS = CI3 / (DP3 / (A * KC3) + 1.0)
                DCPDX = CI3 / A
                TPS = T - (UNBAR('H4TBL', T) * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
                TPS = max(TPS, 1.0)
                continue

            X0 = X0_est
            RIESUM = TRAPP(0.0, 1.0, NPART, 0.0, CPS, GAMMA, BETA, KO)
            print(f"   RIESUM = {RIESUM:.5e}")
            CPSP = CPS
            CMCPO = CMCPN
            CPS = CI3 - A * RIESUM / KC3
            CMCPN = CI3 - CPS
            print(f"   Updated CPS = {CPS:.5e}")

            if LTFLG == 1 and CPS < 0:
                CPS = 0
                continue
            elif CPS < 0.25 * CI3:
                LTFLG = 1
                continue
            if CPS < 0:
                CPS = 0
                continue

            TMTPO = TMTPN
            TMTPN = T - TPS
            print(f"   TPS = {TPS:.2f}, T = {T:.2f}")
            print(f"   TMTPN = {TMTPN:.5e}, TMTPO = {TMTPO:.5e}")
            print(f"   CMCPN = {CMCPN:.5e}, CMCPO = {CMCPO:.5e}")

            if is_converged(TMTPO, TMTPN, CMCPO, CMCPN):
                print(f"\n✅ SATISFACTORY X0 FOUND AFTER {LP1} TRIES")
                print(f"   Final X0 = {X0:.5e}")
                GRAD = DCPDX * DP3
                TGRAD = HC * (T - TPS)
                return GRAD, TGRAD

            if temperature_fluctuates(TPS, TPSP, TPSPP):
                TPS = (TPSP + TPSPP) / 2.0
                H3 = UNBAR('H3TBL', TPS)
                DP3 = DP3F(TPS, D03, P)
                DCPDX = (HC * (T - TPS) - UNBAR('H4TBL', T) * KC4 * CI4) / (H3 * DP3)
                CPS = CI3 - DP3 / KC3 * DCPDX
                CPS = max(CPS, 0)
                continue

            TPSPP = TPSP
            TPSP = TPS
            CPS = 0.2 * CPS + 0.8 * CPSP

        WAF1 += 0.05
        WAF2 = 1.0 - WAF1

    if X0 > 1e-4:
        print("\n⚠️ Proceeding without formal convergence — best effort solution.")
        GRAD = DCPDX * DP3
        TGRAD = HC * (T - TPS)
        return GRAD, TGRAD

    raise RuntimeError(f"UNABLE TO FIND SUITABLE X0 AFTER FOUR TRIES OF 25\nX0 = {X0:.5e}")



def finalize_sgrad_loop(T, P, A, CI1, CI3, CI4, KC3, KC4, KP, DP3, H3, H4, HC,
                        ALPHA2, EN3, X0, CPX, CPOX, PCPOX, RHET, CPS, 
                        eval1, eval2, D03, H3TBL, UNBAR, DP3F, PRINT=0):
    NX = 24
    NX1 = NX + 1
    NXM1 = NX - 1
    DELXOA = (1.0 - X0 / A) / NX
    VNU = -KC3 / DP3
    CTRM = (A * VNU + 1.0) / (A * VNU)

    LP2 = 1
    TMTPN = T - 1.0  # initialize for convergence check
    TMTPO = T
    CMCPN = CI3
    CMCPO = 0.0

    while LP2 <= 50:
        # Generate RHET profile
        XA = X0 / A
        for i in range(NX1):
            if LP2 == 1:
                CPX[i] = (XA - X0 / A) / (1.0 - X0 / A) * CPS
            beta_term = (1.0 - CPX[i] / CI3)
            BETA = -CPX[NX1 - 1] * H3 * DP3 / (KP * T)
            GAMMA = GG / T
            K0 = ALPHA2 * np.exp(-GAMMA) * CI1 ** EN3
            RHET[i] = K0 * CI3 ** (1 - EN3) * CPX[i] ** EN3 * np.exp(
                GAMMA * BETA * beta_term / (1.0 + BETA * beta_term)
            )
            XA += DELXOA

        # Midpoint average
        for i in range(NX):
            CPX[i] = 0.5 * (CPX[i] + CPX[i + 1])
            RHET[i] = 0.5 * (RHET[i] + RHET[i + 1])

        # Solve integral equation at X0
        RR1 = 0.0
        DXL = X0 / A
        DXU = DXL + DELXOA
        for i in range(NX):
            RR1 += RHET[i] * (eval2(DXL, DXU) - CTRM * eval1(DXL, DXU))
            DXL = DXU
            DXU += DELXOA
        CPOX[0] = CI3 - A * A / DP3 * RR1
        if CPOX[0] < 0.0:
            CPOX[0] = 0.0

        # General case (1 < i < NX)
        XOA = X0 / A
        XA = XOA + DELXOA
        INT1 = 1
        for k in range(1, NX):
            R1 = 0.0
            for i in range(INT1):
                R1 += RHET[i] * eval1(XOA, XA)
                XOA = XA
                XA += DELXOA
            R1 *= (1.0 / XOA - CTRM)

            XAD = XA
            XA -= DELXOA
            PS1 = 0.0
            PS2 = 0.0
            for i in range(INT1, NXM1):
                PS1 += RHET[i + 1] * eval2(XA, XAD)
                PS2 += RHET[i + 1] * eval1(XA, XAD)
                XA = XAD
                XAD += DELXOA
            R2 = PS1 - CTRM * PS2
            INT1 += 1
            CPOX[k] = CI3 - A * A / DP3 * (R1 + R2)
            if CPOX[k] < 0.0:
                CPOX[k] = 0.0

        # Boundary at x = A
        RR2 = 0.0
        DXL = X0 / A
        DXU = DXL + DELXOA
        for i in range(NX):
            RR2 += RHET[i] * eval1(DXL, DXU)
            DXL = DXU
            DXU += DELXOA
        CPOX[NX1 - 1] = CI3 - A * A / DP3 * (1.0 - CTRM) * RR2
        if CPOX[NX1 - 1] < 0.0:
            CPOX[NX1 - 1] = 0.0

        # Calculate new temperature TPS
        DCPDX = KC3 / DP3 * (CI3 - CPOX[NX1 - 1])
        TPS = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
        if TPS < 0.0:
            TPS = 1.0

        # Update H3 and DP3 using new TPS
        H3 = UNBAR(H3TBL, 1, TPS, 0.0)
        DP3_new = DP3F(TPS, D03, P)

        # Convergence check after 2 iterations
        if LP2 > 1:
            T_DELTA = abs(TMTPO - TMTPN) / TMTPN
            C_DELTA = abs(CMCPO - CMCPN) / CMCPN if CMCPN != 0 else 1.0
            if T_DELTA <= 0.05 and C_DELTA <= 0.05:
                break

        # Update for next iteration
        TMTPO = TMTPN
        TMTPN = T - TPS
        CMCPO = CMCPN
        CMCPN = CI3 - CPOX[NX1 - 1]

        # Smooth CPX
        for i in range(NX1):
            if LP2 % 5 == 0:
                CPX[i] = 0.5 * (CPOX[i] + PCPOX[i])
            else:
                CPX[i] = 0.8 * CPX[i] + 0.2 * CPOX[i]
            PCPOX[i] = CPOX[i]

        # Update for next iteration
        DP3 = DP3_new
        LP2 += 1

    if LP2 > 50:
        print(f"⚠️ UNABLE TO CONVERGE ON CPS IN 50 TRIES ... CP(X/A) = {CPOX[NX1 - 1]:.5e}")

    GRAD = KC3 * (CI3 - CPOX[NX1 - 1])
    TGRAD = HC * (T - TPS)

    if PRINT == 1:
        print(f"\nCONCENTRATION GRADIENT FOUND AFTER {LP2} TRIES")
        print(f"CP(X) AT PARTICLE SURFACE = {CPOX[NX1 - 1]:.5e}")
        print(f"KC3*(CI3-CPS) = {GRAD:.5e}")
        print(f"HC*(T-TPS)    = {TGRAD:.5e}")

    return GRAD, TGRAD



