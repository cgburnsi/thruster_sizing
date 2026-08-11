"""
Fresh Python translation of Kesten's Fortran SGRAD subroutine.

Translates SGRAD.f + TRAP2.f directly. Fortran line numbers referenced
in comments (e.g. # F:490).

Unit system: Rankine, psia, lb/ft^3, lb/ft^2-s, BTU/...
"""

import numpy as np


# ---------------------------------------------------------------------------
# Inline Fortran functions
# ---------------------------------------------------------------------------

def _dp3f(T, D0, P):
    """DP3F(X,Y,Z) = 14.7*Y/Z*(X/492.)**1.823*(1-EXP(-.0672*Z*492./(14.7*X)))  F:160"""
    T  = max(float(T),  1e-30)
    P  = max(float(P),  1e-30)
    D0 = max(float(D0), 1e-30)
    return 14.7 * D0 / P * (T / 492.0) ** 1.823 * (
        1.0 - np.exp(-0.0672 * P * 492.0 / (14.7 * T))
    )


def _kcf(G, rho, mu, diff, AP):
    """KCF(A,B,C,D,E) = .61*A/B*(C/(B*D))**-.667*(A/(E*C))**-.41  F:180"""
    rho  = max(float(rho),  1e-30)
    mu   = max(float(mu),   1e-30)
    diff = max(float(diff), 1e-30)
    AP   = max(float(AP),   1e-30)
    return 0.61 * G / rho * (mu / (rho * diff)) ** (-0.667) * (G / (AP * mu)) ** (-0.41)


def _eval1(a, b):
    """EVAL1(A,B) = B**3/3. - A**3/3.  F:190"""
    return float(b) ** 3 / 3.0 - float(a) ** 3 / 3.0


def _eval2(a, b):
    """EVAL2(A,B) = B**2/2. - A**2/2.  F:200"""
    return float(b) ** 2 / 2.0 - float(a) ** 2 / 2.0


def _rhetf(ci3, cpx, gamma, beta, k0, order_n):
    """
    RHETF from TRAP2.f F:30
    RHETF(A,B,C,D,E,N) = E*A**(1-N)*B**N*EXP(C*D*(1-B/A)/(1+D*(1-B/A)))
    (N is hardcoded 1 in Fortran TRAPP call; order_n generalises it)
    """
    if ci3 <= 0.0 or cpx <= 0.0 or k0 <= 0.0:
        return 0.0
    theta = 1.0 - cpx / ci3
    denom = max(1.0 + beta * theta, 1e-30)
    return k0 * ci3 ** (1.0 - order_n) * cpx ** order_n * np.exp(
        np.clip(gamma * beta * theta / denom, -200.0, 200.0)
    )


def _cpxf(x, x0a, cps):
    """CPXF(X,Y,Z) = (X-Y)/(1-Y)*Z  F:60 in TRAP2"""
    if x <= x0a:
        return 0.0
    denom = max(1.0 - x0a, 1e-30)
    return (x - x0a) / denom * cps


def _trapp(u, v, npart, x0a, cps, ci3, gamma, beta, k0, order_n):
    """
    Trapezoidal integration of x^2 * rhet(x) from u to v.
    Mirrors TRAPP subroutine (TRAP2.f).
    """
    if v <= u:
        return 0.0
    n    = max(int(npart) - 1, 1)
    part = float(npart)
    h    = (v - u) / part
    uph  = u + h

    cpx1  = _cpxf(u, x0a, cps)
    cpx2  = _cpxf(v, x0a, cps)
    rhet1 = _rhetf(ci3, cpx1, gamma, beta, k0, order_n)
    rhet2 = _rhetf(ci3, cpx2, gamma, beta, k0, order_n)

    # F:160-170 — endpoint half-weights
    trm1 = (u ** 2 * rhet1) / 2.0
    trm2 = (v ** 2 * rhet2) / 2.0   # Fortran has a typo using rhet1 here; rhet2 is correct

    total = 0.0
    for _ in range(n):
        cpx  = _cpxf(uph, x0a, cps)
        rhet = _rhetf(ci3, cpx, gamma, beta, k0, order_n)
        total += uph ** 2 * rhet
        uph += h

    return h * (trm1 + total + trm2)


# ---------------------------------------------------------------------------
# Main SGRAD routine
# ---------------------------------------------------------------------------

def SGRAD(TEMP, G, concentrations, A, AP, target_reaction, target_species,
          cfg, props, pres=None, verbose=False):
    """
    Direct translation of Kesten Fortran SGRAD subroutine.

    Returns (grad, tgrad, diag):
        grad  = KC3 * (CI3 - CPS)  [mass flux, lb/ft^2-s]
        tgrad = HC  * (T  - TPS)   [heat flux, BTU/ft^2-s]
    """
    # --- Only ported for NH3 path ---
    if target_species is None or target_species.name != "NH3":
        return 0.0, 0.0, {"note": "SGRAD only ported for NH3"}

    T = float(TEMP)
    P = float(cfg.PRES if pres is None else pres)

    if T <= 0.0 or P <= 0.0 or A <= 0.0 or AP <= 0.0 or G <= 0.0:
        return 0.0, 0.0, {"note": "degenerate conditions"}

    # --- Species concentrations  F:370 ---
    ci1 = max(float(concentrations.get("H2",   0.0)), 0.0)
    ci2 = max(float(concentrations.get("N2",   0.0)), 0.0)
    ci3 = max(float(concentrations.get("NH3",  0.0)), 0.0)
    ci4 = max(float(concentrations.get("N2H4", 0.0)), 0.0)
    if ci3 <= 0.0:
        return 0.0, 0.0, {"note": "no NH3"}

    rho = max(ci1 + ci2 + ci3 + ci4, 1e-30)    # F:370

    # --- Mixture properties at bulk T  F:380-480 ---
    mu   = max(float(props.VISVST(T)), 1e-30)
    cf1  = props.CFTBL1(T)
    cf2  = props.CFTBL2(T)
    cf3  = props.CFTBL3(T)
    cf4  = props.CFTBL4(T)
    cfbar = (ci1*cf1 + ci2*cf2 + ci3*cf3 + ci4*cf4) / rho   # F:470

    di3 = cfg.DIF3 * 14.7 / P * (T / 492.0) ** 1.823   # F:380 (no correction term)
    di4 = cfg.DIF4 * 14.7 / P * (T / 492.0) ** 1.823   # F:390

    kc3 = _kcf(G, rho, mu, di3, AP)    # F:450
    kc4 = _kcf(G, rho, mu, di4, AP)    # F:460
    hc  = 0.74 * G * cfbar * (G / (AP * mu)) ** (-0.41)   # F:480

    dp3  = _dp3f(T, cfg.DIF3, P)       # F:490
    h4   = props.H4TBL(T)              # F:530
    h3   = props.H3TBL(T)              # F:540
    h3p  = h3
    dp3p = dp3

    # --- EN3 / order for NH3 reaction  F:330 ---
    order_n = float(getattr(cfg, "EN2", 1.0))
    h2_floor = max(ci1, 1e-30)

    # --- Weighted-average loop for X0 (waf1 from 0.80 → 0.95)  F:1390 ---
    WAF1_VALUES = np.arange(0.8, 1.0, 0.05)
    found_x0 = False
    grad  = 0.0
    tgrad = 0.0

    # initialise loop-persistent variables
    tps = T; tpsp = 0.0; tpspp = 0.0
    cps = ci3 / 2.0
    x0 = 0.0; x0p = 0.0
    cmcpn = ci3 - cps
    dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cps)
    gamma = cfg.BGM / max(T, 1e-12)
    beta  = 0.0
    k0    = 0.0

    for waf1 in WAF1_VALUES:
        waf2  = 1.0 - waf1

        # --- Reset per waf attempt  F:1 (LTFLG=0 etc.) ---
        ltflg = 0
        lp1   = 1
        cps   = ci3 / 2.0                          # F:500
        cmcpn = ci3 - cps                          # F:510
        dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cps)  # F:520
        h3    = props.H3TBL(T)
        h3p   = h3
        dp3p  = dp3
        dp3   = _dp3f(T, cfg.DIF3, P)
        tps   = T
        tpsp  = 0.0
        tpspp = 0.0
        x0    = 0.0
        x0p   = 0.0
        tmtpn = 0.0

        while lp1 <= 25:
            # F:560-570 — save last two surface temps
            if lp1 != 1:
                tpspp = tpsp
                tpsp  = tps

            # F:580 — new surface temperature
            tps = T - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
            tps = max(tps, 1.0)                    # F:590

            h3   = props.H3TBL(tps)                # F:600
            dp3  = _dp3f(tps, cfg.DIF3, P)         # F:610
            dp3p = dp3                             # F:620
            h3p  = h3                              # F:630
            tmtpn = T - tps                        # F:640

            # F:650-670 — Gamma, Beta, K0 at surface T
            gamma = cfg.BGM / max(tps, 1e-12)
            beta  = -cps * h3 * dp3 / max(cfg.KP * tps, 1e-30)
            k0    = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3

            # F:690-700 — estimate X0
            x0p = x0
            if dcpdx > 1e-30:
                x0 = A - cps / dcpdx
            else:
                x0 = 0.0

            if x0 < 0.0:                           # F:710-790 — negative X0 branch
                x0  = 0.0
                x0a = 0.0
                denom = max(1.0 + dp3 / max(A * kc3, 1e-30), 1e-30)
                cps   = ci3 / denom
                dcpdx = ci3 / max(A, 1e-30)
                tps   = T - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
                tps   = max(tps, 1.0)
            else:
                x0a = x0 / max(A, 1e-30)

            # F:830 — integrate reaction rate profile
            riesum = _trapp(x0a, 1.0, 50, x0a, cps, ci3, gamma, beta, k0, order_n)

            cpsp  = cps
            cmcpo = cmcpn
            cps   = ci3 - A * riesum / max(kc3, 1e-30)   # F:860

            # F:870-940 — LTFLG branch
            if ltflg != 1:
                if cps <= 0.25 * ci3:              # F:880-890
                    ltflg = 1
                    x00   = waf1 * x0p + waf2 * x0  # F:1320
                    denom = kc3 * A - kc3 * x00
                    if abs(denom) < 1e-30:
                        cps = 0.0
                    else:
                        cps = ci3 / (1.0 + dp3 / denom)  # F:1330
                    dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)  # F:1340
                    cmcpn = ci3 - cps              # F:1350
                    h3    = h3p                    # F:1360
                    lp1  += 1                      # F:1370
                    continue
                cmcpn = ci3 - cps                  # F:950
            else:
                ltflg = 0
                if cps < 0.0:                      # F:920-940
                    cps = 0.0
                    cps = 0.2 * cps + 0.8 * cpsp   # F:1300
                    dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)
                    cmcpn = ci3 - cps
                    h3    = h3p
                    lp1  += 1
                    continue
                cmcpn = ci3 - cps

            # F:960-1000 — update gradient, heat terms, save history
            dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cps)
            tgrad = hc * (T - tps)
            tpspp = tpsp
            tpsp  = tps
            tmtpo = tmtpn

            # F:1020 — recalculate TPS with updated dcpdx
            tps = T - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
            tps = max(tps, 1.0)

            h3    = props.H3TBL(tps)
            dp3   = _dp3f(tps, cfg.DIF3, P)
            tmtpn = T - tps

            gamma = cfg.BGM / max(tps, 1e-12)
            beta  = -cps * h3 * dp3 / max(cfg.KP * tps, 1e-30)
            k0    = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3

            # F:1100-1110 — convergence test
            # Fortran uses signed TMTPN: ABS(TMTPO-TMTPN)/TMTPN <= 0.05
            # When TPS > T, TMTPN = T-TPS < 0 → LHS ≤ 0 ≤ 0.05 → always passes
            temp_ok = tmtpn <= 0.0 or abs(tmtpo - tmtpn) / tmtpn <= 0.05
            conc_ok = abs(cmcpo - cmcpn) / max(abs(cmcpn), 1e-30) <= 0.05
            if temp_ok and conc_ok:
                found_x0 = True
                break

            # F:1120-1130 — temperature oscillation check
            # Fortran: bisect only when TPSP is local extremum (min or max of 3 temps)
            # i.e., temperature is oscillating (TPSP is a turning point)
            if min(tps, tpsp, tpspp) == tpsp or max(tps, tpsp, tpspp) == tpsp:
                # F:1140-1290 — average last two temps (oscillation detected)
                tpspp = tpsp
                tpsp  = tps
                tps   = 0.5 * (tpsp + tpspp)
                h3    = props.H3TBL(tps)
                dp3   = _dp3f(tps, cfg.DIF3, P)
                dp3p  = dp3
                tmtpn = T - tps
                dcpdx = (hc * (T - tps) - h4 * kc4 * ci4) / max(h3 * dp3, 1e-30)
                cpsp  = cps
                cmcpo = cmcpn
                cps   = ci3 - dp3 / max(kc3, 1e-30) * dcpdx
                cps   = max(cps, 0.0)
                cmcpn = ci3 - cps
                lp1  += 1
                continue

            # F:1300 — damp CPS update
            cps   = 0.2 * cps + 0.8 * cpsp
            dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)
            cmcpn = ci3 - cps
            h3    = h3p
            lp1  += 1

        if found_x0:
            break

    if not found_x0:
        grad  = kc3 * max(ci3 - cps, 0.0)
        tgrad = hc  * max(T - tps, 0.0)
        return grad, tgrad, {
            "note": "x0 did not converge",
            "kc3": kc3, "kc4": kc4, "tps": tps, "cps": cps,
        }

    # =======================================================================
    # Phase 2 — integral-equation iteration for concentration profile
    # F:1500
    # =======================================================================
    NX   = 24
    NX1  = NX + 1
    NXM1 = NX - 1

    cpx   = np.zeros(NX1, dtype=float)
    pcpox = np.zeros(NX1, dtype=float)
    cpox  = np.zeros(NX1, dtype=float)
    dx    = np.zeros(NX1, dtype=float)
    rhet  = np.zeros(NX1, dtype=float)

    lp2         = 1
    cpx_surface = cps
    tmtpn       = T - tps
    tmtpo       = tmtpn
    cmcpn       = ci3 - cps
    cmcpo       = cmcpn

    while lp2 <= 300:
        # F:1540
        x0a    = x0 / max(A, 1e-30)
        vnu    = -kc3 / max(dp3, 1e-30)            # F:1550
        ctrm   = (A * vnu + 1.0) / (A * vnu)       # F:1770
        delx0a = (1.0 - x0a) / float(NX)           # F:1620
        xa     = x0a

        # F:1630-1710 — build node values of CPX and RHET
        for i in range(NX1):
            if lp2 == 1:
                cpx[i] = _cpxf(xa, x0a, cps)       # F:1660
            rhet[i] = _rhetf(ci3, cpx[i], gamma, beta, k0, order_n)  # F:1670
            dx[i]   = xa
            xa     += delx0a

        # F:1720-1750 — replace node values with interval midpoints
        for i in range(NX):
            cpx[i]  = 0.5 * (cpx[i]  + cpx[i+1])
            rhet[i] = 0.5 * (rhet[i] + rhet[i+1])

        # ------------------------------------------------------------------
        # F:1780-1860 — CPOX(1): special case at x = x0a
        # ------------------------------------------------------------------
        dxl = x0a
        dxu = dxl + delx0a
        rr1 = 0.0
        for i in range(NX):
            rr1 += rhet[i] * (_eval2(dxl, dxu) - ctrm * _eval1(dxl, dxu))
            dxl  = dxu
            dxu += delx0a
        cpox[0] = max(ci3 - A*A / max(dp3, 1e-30) * rr1, 0.0)    # F:1860

        # ------------------------------------------------------------------
        # F:1880-2130 — CPOX(K) for K = 2 … NX  (k = 1 … NX-1 in Python)
        # Key fix: Fortran DO 773 I=INT1,NXM1 uses RHET(I+1).
        #   In 0-based Python: RHET(I+1) at Fortran-I = rhet[I] at Python-i=I.
        #   Range must be range(int1, NX) to include NXM1=NX-1.
        # ------------------------------------------------------------------
        int1 = 1
        for k in range(1, NX):
            r1  = 0.0;  r2  = 0.0
            ps1 = 0.0;  ps2 = 0.0
            x0a_k = x0 / max(A, 1e-30)
            xa    = x0a_k + delx0a

            # F:1880-1920 — R1 inner integral (INT1 steps)
            for i in range(int1):              # i = 0 … int1-1  ↔  Fortran I=1…INT1
                r1    += rhet[i] * _eval1(x0a_k, xa)
                x0a_k  = xa
                xa    += delx0a
            r1 *= (1.0 / max(x0a_k, 1e-30) - ctrm)   # F:1930

            xad = xa
            xa  = xa - delx0a

            # F:1960-2010 — PS1/PS2 outer integral (INT1 … NXM1 steps)
            # Fortran: DO 773 I=INT1,NXM1  using RHET(I+1)
            # 0-based Python: RHET(I+1)[Fortran] = rhet[I][Python]; range = range(int1, NX)
            for I in range(int1, NX):          # I goes int1 … NX-1 (NXM1) inclusive
                ps1 += rhet[I] * _eval2(xa, xad)
                ps2 += rhet[I] * _eval1(xa, xad)
                xa   = xad
                xad += delx0a

            r2      = ps1 - ctrm * ps2         # F:2020
            cpox[k] = max(ci3 - A*A / max(dp3, 1e-30) * (r1 + r2), 0.0)  # F:2040-2050
            int1   += 1

        # ------------------------------------------------------------------
        # F:2140-2230 — CPOX(NX1): special case at x = A (surface)
        # ------------------------------------------------------------------
        dxl = x0 / max(A, 1e-30)
        dxu = dxl + delx0a
        rr2 = 0.0
        for i in range(NX):
            rr2 += rhet[i] * _eval1(dxl, dxu)
            dxl  = dxu
            dxu += delx0a
        cpox[NX] = max(ci3 - A*A / max(dp3, 1e-30) * (1.0 - ctrm) * rr2, 0.0)  # F:2220

        # F:2240-2310 — update TPS using CPOX at surface
        dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cpox[NX])
        h3p   = h3
        dp3p  = dp3
        tps   = T - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
        tps   = max(tps, 1.0)
        h3    = props.H3TBL(tps)
        dp3   = _dp3f(tps, cfg.DIF3, P)
        tmtpo = tmtpn
        tmtpn = T - tps

        # F:2320-2360 — convergence check (skip first pass)
        if lp2 > 1:
            cmcpo = cmcpn
            cmcpn = ci3 - cpox[NX]
            # F:2350: ABS(TMTPO-TMTPN)/TMTPN ≤ 0.05 — signed denominator
            temp_ok = tmtpn <= 0.0 or abs(tmtpo - tmtpn) / tmtpn <= 0.05
            conc_ok = abs(cmcpo - cmcpn) / max(abs(cmcpn), 1e-30) <= 0.05
            if temp_ok and conc_ok:
                cpx_surface = cpox[NX]
                break

        # F:2370-2530 — update CPX profile for next pass
        for i in range(NX1):
            if lp2 % 5 != 0:                   # F:2380 — weighted average
                cpx[i] = 0.8 * cpx[i] + 0.2 * cpox[i]   # F:2390
            else:                               # F:2410 — smooth every 5th pass
                cpx[i] = 0.5 * (cpox[i] + pcpox[i])
            pcpox[i] = cpox[i]                 # F:2420

        cmcpn = ci3 - cpx[NX]
        dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cpx[NX])
        tps   = T - (h4 * kc4 * ci4 + h3p * dp3p * dcpdx) / max(hc, 1e-30)
        tps   = max(tps, 1.0)
        h3    = props.H3TBL(tps)
        dp3   = _dp3f(tps, cfg.DIF3, P)
        tmtpn = T - tps

        gamma = cfg.BGM / max(tps, 1e-12)              # F:2600
        beta  = -cpx[NX] * h3 * dp3 / max(cfg.KP * tps, 1e-30)   # F:2610
        k0    = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3  # F:2620

        cpx_surface = cpx[NX]
        lp2 += 1

    # F:2710-2720 — final output
    grad  = kc3 / max(dp3, 1e-30) * (ci3 - cpx_surface) * dp3   # = kc3*(ci3-cps)
    tgrad = hc * (T - tps)

    diag = {
        "tps":   tps,
        "cps":   cpx_surface,
        "x0":    x0,
        "x0a":   x0 / max(A, 1e-30),
        "kc3":   kc3,
        "kc4":   kc4,
        "hc":    hc,
        "gamma": gamma,
        "beta":  beta,
        "k0":    k0,
        "lp2":   lp2,
    }
    if verbose:
        print(f"[SGRAD] x0={x0:.4e}  cps={cpx_surface:.4e}  tps={tps:.2f}  lp2={lp2}")
    return grad, tgrad, diag
