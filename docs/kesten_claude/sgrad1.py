import numpy as np


# ---------------------------------------------------------------------------
# Inline Fortran statement functions
# ---------------------------------------------------------------------------

def _dp3f(T, D0, P):
    """DP3F(X,Y,Z) = 14.7*Y/Z*(X/492.)**1.823*(1-EXP(-.0672*Z*492./(14.7*X)))"""
    T  = max(float(T),  1e-30)
    P  = max(float(P),  1e-30)
    D0 = max(float(D0), 1e-30)
    return 14.7 * D0 / P * (T / 492.0)**1.823 * (
        1.0 - np.exp(-0.0672 * P * 492.0 / (14.7 * T))
    )

def _kcf(G, rho, mu, diff, AP):
    """KCF(A,B,C,D,E) = .61*A/B*(C/(B*D))**-.667*(A/(E*C))**-.41"""
    rho  = max(float(rho),  1e-30)
    mu   = max(float(mu),   1e-30)
    diff = max(float(diff), 1e-30)
    AP   = max(float(AP),   1e-30)
    return 0.61 * G / rho * (mu / (rho * diff))**(-0.667) * (G / (AP * mu))**(-0.41)

def _eval1(a, b):
    """EVAL1(A,B) = B**3/3. - A**3/3."""
    return float(b)**3 / 3.0 - float(a)**3 / 3.0

def _eval2(a, b):
    """EVAL2(A,B) = B**2/2. - A**2/2."""
    return float(b)**2 / 2.0 - float(a)**2 / 2.0

def _cpxf(x, x0a, cps):
    """CPXF(X,Y,Z) = (X-Y)/(1-Y)*Z"""
    if x <= x0a:
        return 0.0
    return (x - x0a) / max(1.0 - x0a, 1e-30) * cps

def _rhetf(ci3, cpx, gamma, beta, k0, N):
    """RHETF: reaction rate at concentration cpx"""
    if ci3 <= 0.0 or cpx <= 0.0 or k0 <= 0.0:
        return 0.0
    theta = 1.0 - cpx / ci3
    denom = max(1.0 + beta * theta, 1e-30)
    return k0 * ci3**(1.0 - N) * cpx**N * np.exp(
        np.clip(gamma * beta * theta / denom, -200.0, 200.0)
    )

def _trapp(u, v, npart, x0a, cps, ci3, gamma, beta, k0, N):
    """Trapezoidal integration of x^2 * rhet(x) from u to v (TRAP2.f)."""
    if v <= u:
        return 0.0
    n   = max(int(npart) - 1, 1)
    h   = (v - u) / float(npart)
    uph = u + h

    cpx1  = _cpxf(u, x0a, cps)
    cpx2  = _cpxf(v, x0a, cps)
    rhet1 = _rhetf(ci3, cpx1, gamma, beta, k0, N)
    rhet2 = _rhetf(ci3, cpx2, gamma, beta, k0, N)
    trm1  = u**2 * rhet1 / 2.0
    trm2  = v**2 * rhet2 / 2.0

    total = 0.0
    for _ in range(n):
        cpx   = _cpxf(uph, x0a, cps)
        rhet  = _rhetf(ci3, cpx, gamma, beta, k0, N)
        total += uph**2 * rhet
        uph   += h

    return h * (trm1 + total + trm2)


# ---------------------------------------------------------------------------
# Main SGRAD routine
# ---------------------------------------------------------------------------

def SGRAD(TEMP, G, concentrations, A, AP, target_reaction, target_species,
          cfg, props, pres=None, verbose=False):
    """
    Direct port of Kesten SGRAD Fortran subroutine.

    Returns (GRAD, TGRAD, diag):
        GRAD  = KC3 * (CI3 - CPS)  [mass flux, lb/ft^2-s]
        TGRAD = HC  * (T  - TPS)   [heat flux, BTU/ft^2-s]
    """
    if target_species is None or target_species.name != "NH3":
        return 0.0, 0.0, {"note": "SGRAD only ported for NH3"}

    T = float(TEMP)
    P = float(cfg.PRES if pres is None else pres)

    if T <= 0.0 or P <= 0.0 or A <= 0.0 or AP <= 0.0 or G <= 0.0:
        return 0.0, 0.0, {"note": "degenerate conditions"}

    # -------------------------------------------------------------------------
    # INITIALIZATION
    # -------------------------------------------------------------------------
    WAF1, WAF2 = 0.8, 0.2
    NPART = 50
    NX    = 24
    NX1   = NX + 1

    # Species concentrations
    CI1 = max(float(concentrations.get("H2",   0.0)), 0.0)
    CI2 = max(float(concentrations.get("N2",   0.0)), 0.0)
    CI3 = max(float(concentrations.get("NH3",  0.0)), 0.0)
    CI4 = max(float(concentrations.get("N2H4", 0.0)), 0.0)
    if CI3 <= 0.0:
        return 0.0, 0.0, {"note": "no NH3"}

    D03, D04 = cfg.DIF3, cfg.DIF4
    N        = float(getattr(cfg, "EN1", 1.0))   # reaction order (=1.0)
    EN3      = cfg.EN3                            # H2 inhibition exponent
    h2_floor = max(CI1, 1e-30)

    RHO = max(CI1 + CI2 + CI3 + CI4, 1e-30)
    DI3 = D03 * 14.7 / P * (T / 492.0)**1.823
    DI4 = D04 * 14.7 / P * (T / 492.0)**1.823

    MU   = max(float(props.VISVST(T)), 1e-30)
    CF4  = props.CFTBL4(T)
    CF3  = props.CFTBL3(T)
    CF2  = props.CFTBL2(T)
    CF1  = props.CFTBL1(T)

    KC3  = _kcf(G, RHO, MU, DI3, AP)
    KC4  = _kcf(G, RHO, MU, DI4, AP)

    CFBAR = (CI1*CF1 + CI2*CF2 + CI3*CF3 + CI4*CF4) / RHO
    HC    = 0.74 * G * CFBAR * (G / (AP * MU))**(-0.41)

    DP3 = _dp3f(T, D03, P)

    # Initial guesses
    CPS   = CI3 / 2.0
    CMCPN = CI3 - CPS
    DCPDX = KC3 / max(DP3, 1e-30) * (CI3 - CPS)

    H4 = props.H4TBL(T)
    H3 = props.H3TBL(T)

    X0  = 0.0
    X0A = 0.0
    GRAD  = 0.0
    TGRAD = 0.0

    # -------------------------------------------------------------------------
    # OUTER WAF RETRY LOOP  (label 1 in Fortran)
    # Retries up to 4 times with increasing WAF1 weight (0.80 → 0.95)
    # -------------------------------------------------------------------------
    while WAF1 <= 0.95:

        LP1   = 0
        LTFLG = 0
        TPSP  = 0.0
        TPSPP = 0.0
        X0    = 0.0
        X0A   = 0.0
        TMTPN = 0.0

        # Reset to initial guesses for each WAF attempt
        CPS   = CI3 / 2.0
        CMCPN = CI3 - CPS
        DCPDX = KC3 / max(DP3, 1e-30) * (CI3 - CPS)
        H3    = props.H3TBL(T)
        DP3   = _dp3f(T, D03, P)

        # Compute initial TPS guess
        TPS = T - (H4*KC4*CI4 + H3*DP3*DCPDX) / max(HC, 1e-30)
        TPS = max(TPS, 1.0)
        H3  = props.H3TBL(TPS)
        DP3 = _dp3f(TPS, D03, P)
        DP3P, H3P = DP3, H3

        # -----------------------------------------------------------------
        # PHASE 1 INNER LOOP  — find X0 and CPS  (label 40/42 in Fortran)
        # Up to 25 iterations
        # -----------------------------------------------------------------
        while LP1 < 25:

            GAMMA = cfg.BGM / max(TPS, 1e-12)
            BETA  = -CPS * H3 * DP3 / max(cfg.KP * TPS, 1e-30)
            K0    = cfg.ALPHA2 * np.exp(-GAMMA) * h2_floor**EN3

            # Linear extrapolation for X0
            X0P = X0
            if DCPDX > 1e-30:
                X0 = A - CPS / DCPDX
            else:
                X0 = 0.0

            if X0 < 0:
                # Negative X0: force to zero and recalculate (label 11)
                X0, X0A = 0.0, 0.0
                CPS   = CI3 / (max(DP3 / max(A * KC3, 1e-30), 0.0) + 1.0)
                DCPDX = CI3 / max(A, 1e-30)
                TPS   = T - (H4*KC4*CI4 + H3*DP3*DCPDX) / max(HC, 1e-30)
                TPS   = max(TPS, 1.0)
                if verbose:
                    print("WARNING: negative X0 during iteration")
            else:
                X0A = X0 / max(A, 1e-30)

            # Integrate for CPS using trapezoidal rule (label 12 / TRAPP)
            RIESUM = _trapp(X0A, 1.0, NPART, X0A, CPS, CI3, GAMMA, BETA, K0, N)

            CPSP  = CPS
            CMCPO = CMCPN
            CPS   = CI3 - A * RIESUM / max(KC3, 1e-30)

            # ------------------------------------------------------------------
            # LTFLG branch (Fortran labels 80/84/130/13)
            # ------------------------------------------------------------------
            if LTFLG == 0 and CPS < 0.25 * CI3:
                # Label 81 → 22: weighted average X0 path, then label 53 → 42
                LTFLG = 1
                X00   = WAF1*X0P + WAF2*X0
                denom = KC3*A - KC3*X00
                CPS   = CI3 / (1.0 + DP3 / max(denom, 1e-30)) if abs(denom) > 1e-30 else 0.0
                DCPDX = KC3 / max(DP3P, 1e-30) * (CI3 - CPS)
                CMCPN = CI3 - CPS
                H3    = H3P
                LP1  += 1
                continue

            elif LTFLG == 1 and CPS < 0.0:
                # Label 89 → 46 → 53 → 42: CPS=0 damped path
                LTFLG = 0
                CPS   = 0.0
                CPS   = 0.2*CPS + 0.8*CPSP
                DCPDX = KC3 / max(DP3P, 1e-30) * (CI3 - CPS)
                CMCPN = CI3 - CPS
                H3    = H3P
                LP1  += 1
                continue

            else:
                # Label 130 → 13: normal convergence path
                # Reached when: (LTFLG==0 AND CPS >= 0.25*CI3) OR (LTFLG==1 AND CPS >= 0)
                if LTFLG == 1:
                    LTFLG = 0
                CMCPN = CI3 - CPS

                DCPDX = KC3 / max(DP3, 1e-30) * (CI3 - CPS)
                GRAD  = DCPDX * DP3
                TGRAD = HC * (T - TPS)
                TPSPP = TPSP
                TPSP  = TPS
                TMTPO = TMTPN

                # Recalculate TPS (label 51)
                TPS   = T - (H4*KC4*CI4 + H3*DP3*DCPDX) / max(HC, 1e-30)
                TPS   = max(TPS, 1.0)
                H3    = props.H3TBL(TPS)
                DP3   = _dp3f(TPS, D03, P)
                TMTPN = T - TPS

                GAMMA = cfg.BGM / max(TPS, 1e-12)
                BETA  = -CPS * H3 * DP3 / max(cfg.KP * TPS, 1e-30)
                K0    = cfg.ALPHA2 * np.exp(-GAMMA) * h2_floor**EN3

                # --- Convergence test (labels 1100-1110) ---
                # Fortran uses signed TMTPN: when TPS > T, TMTPN < 0 → always passes
                temp_ok = TMTPN <= 0.0 or abs(TMTPO - TMTPN) / TMTPN <= 0.05
                conc_ok = abs(CMCPO - CMCPN) / max(abs(CMCPN), 1e-30) <= 0.05
                if temp_ok and conc_ok:
                    break  # Phase 1 converged → go to Phase 2

                # --- Oscillation detection (labels 43, 60, 71) ---
                # Bisect only when TPSP is local extremum (oscillating)
                if (min(TPS, TPSP, TPSPP) == TPSP or
                        max(TPS, TPSP, TPSPP) == TPSP):
                    TPSPP = TPSP
                    TPSP  = TPS
                    TPS   = (TPSP + TPSPP) / 2.0
                    H3    = props.H3TBL(TPS)
                    DP3   = _dp3f(TPS, D03, P)
                    DP3P  = DP3
                    TMTPN = T - TPS
                    DCPDX = (HC*(T-TPS) - H4*KC4*CI4) / max(H3*DP3, 1e-30)
                    CPSP  = CPS
                    CMCPO = CMCPN
                    CPS   = CI3 - DP3 / max(KC3, 1e-30) * DCPDX
                    CPS   = max(CPS, 0.0)
                    CMCPN = CI3 - CPS
                else:
                    # Monotone trend — damp CPS (label 46 → 53)
                    CPS   = 0.2*CPS + 0.8*CPSP
                    DCPDX = KC3 / max(DP3P, 1e-30) * (CI3 - CPS)
                    CMCPN = CI3 - CPS
                    H3    = H3P

            LP1 += 1

        else:
            # LP1 hit 25 — try new WAF weights (label 44 → 1390)
            WAF1 += 0.05
            WAF2  = 1.0 - WAF1
            continue

        break  # Phase 1 converged, exit WAF loop

    else:
        # All WAF attempts failed (label 99)
        GRAD  = KC3 * max(CI3 - CPS, 0.0)
        TGRAD = HC  * max(T - TPS, 0.0)
        return GRAD, TGRAD, {
            "note": "x0 did not converge",
            "kc3": KC3, "kc4": KC4, "tps": TPS, "cps": CPS,
        }

    if verbose:
        print(f"Satisfactory X0 found after {LP1} tries.  X0 = {X0:.5e}")

    # =========================================================================
    # PHASE 2 — compute full CP(X/A) spatial profile  (label 131 in Fortran)
    # =========================================================================
    LP2   = 1
    CPOX  = [0.0] * (NX1 + 1)
    CPX   = [0.0] * (NX1 + 1)
    PCPOX = [0.0] * (NX1 + 1)
    DX    = [0.0] * NX1

    TMTPN = T - TPS
    TMTPO = TMTPN
    CMCPN = CI3 - CPS
    CMCPO = CMCPN

    while LP2 <= 50:

        # --- label 291: top of LP2 loop ---
        X0A    = X0 / max(A, 1e-30)
        VNU    = -KC3 / max(DP3, 1e-30)          # VNU < 0 always
        AV     = A * VNU                          # AV < 0 always; keep sign
        CTRM   = (AV + 1.0) / AV                 # = 1 - DP3/(A*KC3)
        DELX0A = (1.0 - X0A) / float(NX)

        # ------------------------------------------------------------------
        # DO 770: build CPX and RHET profiles at NX+1 points
        # ------------------------------------------------------------------
        RHET_pt = [0.0] * NX1
        XA = X0A
        for i in range(NX1):
            if LP2 == 1:
                CPX[i] = _cpxf(XA, X0A, CPS)
            RHET_pt[i] = _rhetf(CI3, CPX[i], GAMMA, BETA, K0, N)
            DX[i] = XA
            XA   += DELX0A

        # ------------------------------------------------------------------
        # DO 771: replace point values with interval midpoint averages
        # ------------------------------------------------------------------
        for i in range(NX):
            CPX[i]     = (CPX[i]     + CPX[i+1])     / 2.0
            RHET_pt[i] = (RHET_pt[i] + RHET_pt[i+1]) / 2.0

        # ------------------------------------------------------------------
        # DO 377: CPOX[0] — special case at X = X0
        # ------------------------------------------------------------------
        DXL, DXU = X0A, X0A + DELX0A
        RR1 = 0.0
        for i in range(NX):
            RR1 += RHET_pt[i] * (_eval2(DXL, DXU) - CTRM * _eval1(DXL, DXU))
            DXL  = DXU
            DXU += DELX0A
        CPOX[0] = max(CI3 - A*A / max(DP3, 1e-30) * RR1, 0.0)

        # ------------------------------------------------------------------
        # Loop 769: CPOX[1] .. CPOX[NX-1] — general interior points
        # Fortran: K=2..NX, INT1=1..NX-1
        # 0-based: CPOX[K-1] = CPOX[1..NX-1]
        # ------------------------------------------------------------------
        INT1      = 1
        K         = 2
        X0A_k     = X0 / max(A, 1e-30)
        XA        = X0A_k + DELX0A
        R1 = R2 = PS1 = PS2 = 0.0

        while K <= NX:
            # DO 772: accumulate left integral R1 (INT1 steps)
            for i in range(INT1):
                R1    += RHET_pt[i] * _eval1(X0A_k, XA)
                X0A_k  = XA
                XA    += DELX0A
            R1 *= (1.0 / max(X0A_k, 1e-30) - CTRM)

            XAD = XA
            XA -= DELX0A

            # DO 773: accumulate right integrals PS1, PS2
            # Fortran: DO 773 I=INT1,NXM1  using RHET(I+1)
            # 0-based: RHET(I+1)[1-based] = RHET_pt[I][0-based]; range = range(INT1, NX)
            for I in range(INT1, NX):
                PS1 += RHET_pt[I] * _eval2(XA, XAD)
                PS2 += RHET_pt[I] * _eval1(XA, XAD)
                XA   = XAD
                XAD += DELX0A

            R2 = PS1 - CTRM * PS2

            INT1      += 1
            CPOX[K-1]  = max(CI3 - A*A / max(DP3, 1e-30) * (R1 + R2), 0.0)

            # Reset accumulators for next K
            X0A_k = X0 / max(A, 1e-30)
            XA    = X0A_k + DELX0A
            R1 = R2 = PS1 = PS2 = 0.0
            K += 1

        # ------------------------------------------------------------------
        # DO 378: CPOX[NX] — special case at X = A
        # ------------------------------------------------------------------
        DXL = X0 / max(A, 1e-30)
        DXU = DXL + DELX0A
        RR2 = 0.0
        for i in range(NX):
            RR2 += RHET_pt[i] * _eval1(DXL, DXU)
            DXL  = DXU
            DXU += DELX0A
        CPOX[NX] = max(CI3 - A*A / max(DP3, 1e-30) * (1.0 - CTRM) * RR2, 0.0)

        # ------------------------------------------------------------------
        # Update TPS from new surface concentration
        # ------------------------------------------------------------------
        DCPDX = KC3 / max(DP3, 1e-30) * (CI3 - CPOX[NX])
        H3P, DP3P = H3, DP3
        TPS   = T - (H4*KC4*CI4 + H3*DP3*DCPDX) / max(HC, 1e-30)
        TPS   = max(TPS, 1.0)
        H3    = props.H3TBL(TPS)
        DP3   = _dp3f(TPS, D03, P)
        TMTPO = TMTPN
        TMTPN = T - TPS

        # ------------------------------------------------------------------
        # Convergence check (skip on first pass — label 33)
        # ------------------------------------------------------------------
        if LP2 > 1:
            CMCPO = CMCPN
            CMCPN = CI3 - CPOX[NX]
            temp_ok = TMTPN <= 0.0 or abs(TMTPO - TMTPN) / TMTPN <= 0.05
            conc_ok = abs(CMCPO - CMCPN) / max(abs(CMCPN), 1e-30) <= 0.05
            if temp_ok and conc_ok:
                CPS = CPOX[NX]
                break

        # ------------------------------------------------------------------
        # DO 55: update CPX profile for next pass
        # ------------------------------------------------------------------
        for i in range(NX1):
            if LP2 % 5 != 0:
                CPX[i] = 0.8*CPX[i] + 0.2*CPOX[i]
            else:
                CPX[i] = (CPOX[i] + PCPOX[i]) / 2.0
            PCPOX[i] = CPOX[i]

        CMCPN = CI3 - CPX[NX]
        DCPDX = KC3 / max(DP3P, 1e-30) * (CI3 - CPX[NX])
        TPS   = T - (H4*KC4*CI4 + H3P*DP3P*DCPDX) / max(HC, 1e-30)
        TPS   = max(TPS, 1.0)
        H3    = props.H3TBL(TPS)
        DP3   = _dp3f(TPS, D03, P)
        TMTPO = TMTPN
        TMTPN = T - TPS

        GAMMA = cfg.BGM / max(TPS, 1e-12)
        BETA  = -CPX[NX] * H3 * DP3 / max(cfg.KP * TPS, 1e-30)
        K0    = cfg.ALPHA2 * np.exp(-GAMMA) * h2_floor**EN3
        CPS   = CPX[NX]
        LP2  += 1

    else:
        if verbose:
            print(f"Unable to converge on CPS in 50 tries. CP(X/A) = {CPOX[NX]:.5e}")

    # -------------------------------------------------------------------------
    # FINAL OUTPUT  (label 88)
    # -------------------------------------------------------------------------
    GRAD  = DCPDX * DP3
    TGRAD = HC * (T - TPS)

    diag = {
        "tps":   TPS,
        "cps":   CPS,
        "x0":    X0,
        "x0a":   X0 / max(A, 1e-30),
        "kc3":   KC3,
        "kc4":   KC4,
        "hc":    HC,
        "gamma": GAMMA,
        "beta":  BETA,
        "k0":    K0,
        "lp1":   LP1,
        "lp2":   LP2,
    }

    if verbose:
        print(f"[SGRAD] x0={X0:.4e}  cps={CPS:.4e}  tps={TPS:.2f}  lp1={LP1}  lp2={LP2}")

    return GRAD, TGRAD, diag
