# SGRAD.f — Codex Translation Briefing

## What this code is and what you need to produce

`SGRAD` is a Fortran IV subroutine that computes two scalar outputs:

- **`GRAD`** — the mass transfer flux at a reacting particle surface: `KC3 * (CI3 - CPS)`
- **`TGRAD`** — the heat flux at the same surface: `HC * (T - TPS)`

It does this by iteratively solving for two unknowns that are coupled to each other:
- **`CPS`** — the concentration of species 3 at the particle surface
- **`TPS`** — the temperature at the particle surface

The coupling is the core difficulty: CPS depends on TPS through the reaction rate, and TPS depends on CPS through the heat balance. Neither can be solved without the other, so the code iterates until both are self-consistent.

Your task is to translate this into clean, working Python. The pseudocode below is a faithful
structural map of the original. Follow it exactly — do not simplify the loop logic or convergence
tests, as they encode hard-won numerical stability that is not obvious from first principles.

---

## Global inputs (these arrive via COMMON blocks in Fortran — treat as function arguments or a config object)

```
# Physical properties
P       = PRES          # system pressure (psia)
T       = TEMP          # bulk gas temperature (°R)
A       = A             # particle radius
G       = G             # mass flux (lb/ft²·s)
AP      = AP            # particle surface area

# Species concentrations (mass fractions or molar)
CI1, CI2, CI3, CI4      # bulk concentrations of species 1-4
                        # CI3 is the key reacting species

# Diffusivities at reference conditions
DIF3, DIF4              # binary diffusivities for species 3 and 4

# Kinetic and transport parameters
BGM                     # activation temperature (gamma numerator)
KP                      # thermal conductivity of particle
ALPH2                   # pre-exponential factor
EN3                     # reaction order

# Heat transfer params
H4                      # enthalpy of species 4 (looked up from table at T)
KC4                     # mass transfer coefficient for species 4

# Lookup tables (1D interpolation arrays)
VISVST[30]              # viscosity vs temperature
CFTBL1..4[34]           # heat capacity vs temperature for each species
H4TBL[40], H3TBL[40]   # enthalpy vs temperature for species 3 and 4
```

---

## Helper functions — implement these first

```python
def DP3F(temp, D03, pressure):
    """Diffusivity of species 3 corrected for T and P."""
    return 14.7 * D03 / pressure * (temp / 492.0) ** 1.823

def KCF(G, RHO, MU, DI, AP):
    """Mass transfer coefficient."""
    return 0.61 * G / RHO * (MU / (RHO * DI)) ** -0.667 * (G / (AP * MU)) ** -0.41

def EVAL1(a, b):
    """Analytic integral of x^2 from a to b (i.e. b³/3 - a³/3)."""
    return b**3 / 3.0 - a**3 / 3.0

def EVAL2(a, b):
    """Analytic integral of x from a to b (i.e. b²/2 - a²/2)."""
    return b**2 / 2.0 - a**2 / 2.0

def unbar(table, T):
    """
    1D linear interpolation into a lookup table.
    The original UNBAR subroutine takes (table, flag, T, 0., result, KK).
    Implement as a standard linear interpolation — the table is evenly spaced
    in temperature. You will need to know the temperature range and spacing
    for each table. Treat these as known constants from the original code's
    COMMON /FTZ/ block.
    """
    pass  # implement with numpy.interp or equivalent

def TRAPP(x0a, upper, npart, beta, gamma, k0, ci3, n):
    """
    Trapezoidal integration of the reaction rate function from x0a to upper
    over npart intervals. Called in Phase 1 to get RIESUM, which is then
    used to update CPS.

    The integrand involves the linear CP profile and the RHET expression:
        integrand(x) = K0 * CI3^(1-N) * cp(x)^N
                       * exp(GAMMA * BETA * (1 - cp(x)/CI3)
                             / (1 + BETA * (1 - cp(x)/CI3)))
    where cp(x) is the linear profile: cp(x) = (x - x0a) / (1 - x0a) * CPS

    Returns RIESUM.
    """
    pass
```

---

## Top-level structure

The function has three nested levels of iteration. Read this carefully before writing any code.

```
SGRAD
│
├── Initialization
│     Compute RHO, DI3, DI4, MU, CF1-CF4, KC3, KC4, CFBAR, HC, DP3
│     Initial guess: CPS = CI3 / 2
│
├── OUTER RETRY LOOP  (WAF loop)
│     WAF1 starts at 0.80, increments by 0.05 each retry
│     WAF2 = 1 - WAF1
│     Abort with error if WAF1 > 0.95  (more than 4 retries)
│     Resets LP1 = 0 on each retry
│
│   ├── PHASE 1 LOOP  (LP1, max 25 iterations)
│   │     Finds self-consistent CPS, TPS, X0 using a LINEAR cp profile
│   │     Convergence: both temperature and concentration residuals < 5%
│   │     On non-convergence after 25 iters: increment WAF1, retry outer loop
│   │
│   └── [Phase 1 converged] → fall through to Phase 2
│
└── PHASE 2 LOOP  (LP2, max 50 iterations)
      Finds the full spatial CP(X/A) profile using the integral equation
      Updates CPS, TPS, K0, GAMMA, BETA each pass
      Convergence: same 5% test on temperature and concentration
      On non-convergence after 50 iters: print warning, return best estimate
```

---

## Phase 1 — detailed logic

**Goal:** Find X0 (flame front location) and CPS (surface concentration) such that the
heat balance and mass balance are simultaneously satisfied, using a linear approximation
for the CP(X/A) profile.

```python
LP1   = 0
LTFLG = 0       # flag: has CPS gone below 0.25*CI3?
TPSP  = 0.0     # TPS from previous iteration
TPSPP = 0.0     # TPS from two iterations ago

# First TPS estimate
TPS   = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
TPS   = max(TPS, 1.0)
H3    = unbar(H3TBL, TPS)
DP3   = DP3F(TPS, D03, P)
DP3P  = DP3
H3P   = H3
TMTPN = T - TPS

while LP1 < 25:

    GAMMA = BGM / TPS
    BETA  = -CPS * H3 * DP3 / (KP * TPS)
    K0    = ALPH2 * exp(-GAMMA) * CI1**EN3

    # Extrapolate X0
    X0P = X0
    X0  = A - CPS / DCPDX

    if X0 < 0:
        X0, X0A = 0.0, 0.0
        CPS   = CI3 / (DP3 / (A * KC3) + 1.0)
        DCPDX = CI3 / A
        TPS   = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
        TPS   = max(TPS, 1.0)
        print("WARNING: negative X0 during iteration")
        # Fall through to TRAPP call below

    # Integrate (always called regardless of X0 sign path)
    RIESUM = TRAPP(X0A, 1.0, NPART, ...)

    CPSP  = CPS
    CMCPO = CMCPN
    CPS   = CI3 - A * RIESUM / KC3

    # ---- CPS routing logic (LTFLG state machine) ----
    if LTFLG == 0 and CPS < 0.25 * CI3:
        # CPS has gone too low — use weighted average X0 instead
        LTFLG = 1
        X00   = WAF1 * X0P + WAF2 * X0
        CPS   = CI3 / (1.0 + DP3 / (KC3 * A - KC3 * X00))
        DCPDX = KC3 / DP3P * (CI3 - CPS)
        CMCPN = CI3 - CPS
        H3    = H3P
        LP1  += 1
        continue

    elif LTFLG == 1:
        LTFLG = 0
        if CPS <= 0.0:
            CPS   = 0.0
            CPS   = 0.2 * CPS + 0.8 * CPSP   # damped update
            DCPDX = KC3 / DP3P * (CI3 - CPS)
            CMCPN = CI3 - CPS
            H3    = H3P
            LP1  += 1
            continue

    # ---- Normal update path ----
    CMCPN = CI3 - CPS
    DCPDX = KC3 / DP3 * (CI3 - CPS)
    GRAD  = DCPDX * DP3
    TGRAD = HC * (T - TPS)

    TPSPP = TPSP
    TPSP  = TPS
    TMTPO = TMTPN

    TPS   = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
    TPS   = max(TPS, 1.0)
    H3    = unbar(H3TBL, TPS)
    DP3   = DP3F(TPS, D03, P)
    TMTPN = T - TPS
    GAMMA = BGM / TPS
    BETA  = -CPS * H3 * DP3 / (KP * TPS)
    K0    = ALPH2 * exp(-GAMMA) * CI1**EN3

    # ---- Convergence test ----
    if (abs(TMTPO - TMTPN) / TMTPN < 0.05 and
        abs(CMCPO - CMCPN) / CMCPN < 0.05):
        break   # Phase 1 done

    # ---- Oscillation detection ----
    # If TPS is not monotone across last 3 values, average to damp
    if (min(TPS, TPSP, TPSPP) == TPSP or
        max(TPS, TPSP, TPSPP) == TPSP):
        TPS   = (TPSP + TPSPP) / 2.0
        H3    = unbar(H3TBL, TPS)
        DP3   = DP3F(TPS, D03, P)
        DP3P  = DP3
        TMTPN = T - TPS
        DCPDX = (HC * (T - TPS) - H4 * KC4 * CI4) / (H3 * DP3)
        CPSP  = CPS
        CMCPO = CMCPN
        CPS   = CI3 - DP3 / KC3 * DCPDX
        CPS   = max(CPS, 0.0)
        CMCPN = CI3 - CPS

    LP1 += 1

else:
    # LP1 exhausted — trigger WAF retry
    WAF1 += 0.05
    WAF2  = 1.0 - WAF1
    if WAF1 > 0.95:
        raise RuntimeError(f"Cannot find X0 after 4 retries. X0 = {X0:.5e}")
    continue   # restart outer WAF loop
```

---

## Phase 2 — detailed logic

**Goal:** Compute the full spatial concentration profile CPOX(I) over the particle,
from X=X0 (flame front) to X=A (surface), using the integral equation of the
boundary value problem. Iterate until the profile is self-consistent with TPS.

```python
LP2   = 1
NX    = 24
NX1   = NX + 1
CPOX  = [0.0] * NX1    # calculated profile this pass
CPX   = [0.0] * NX1    # working / smoothed profile (starts linear from Phase 1)
PCPOX = [0.0] * NX1    # previous calculated profile (for every-5th-pass smoothing)

while LP2 <= 50:

    X0A    = X0 / A
    VNU    = -KC3 / DP3
    DELX0A = (1.0 - X0A) / NX
    CTRM   = (A * VNU + 1.0) / (A * VNU)

    # ------------------------------------------------------------------
    # Step A: Build CPX and RHET profiles at NX+1 points
    # ------------------------------------------------------------------
    DX      = [0.0] * NX1
    RHET_pt = [0.0] * NX1
    XA = X0A
    for i in range(NX1):
        if LP2 == 1:
            CPX[i] = (XA - X0A) / (1.0 - X0A) * CPS   # linear on first pass
        RHET_pt[i] = (K0 * CI3**(1 - N) * CPX[i]**N
                      * exp(GAMMA * BETA * (1.0 - CPX[i] / CI3)
                            / (1.0 + BETA * (1.0 - CPX[i] / CI3))))
        DX[i] = XA
        XA   += DELX0A

    # Replace point values with interval midpoint averages
    for i in range(NX):
        CPX[i]     = (CPX[i]     + CPX[i + 1])     / 2.0
        RHET_pt[i] = (RHET_pt[i] + RHET_pt[i + 1]) / 2.0

    # ------------------------------------------------------------------
    # Step B: CPOX[0] — special case at X = X0 (left boundary)
    # This uses only the left-sum integral (DO 377 in Fortran)
    # ------------------------------------------------------------------
    DXL, DXU = X0A, X0A + DELX0A
    RR1 = 0.0
    for i in range(NX):
        RR1 += RHET_pt[i] * (EVAL2(DXL, DXU) - CTRM * EVAL1(DXL, DXU))
        DXL  = DXU
        DXU += DELX0A
    CPOX[0] = max(CI3 - A * A / DP3 * RR1, 0.0)

    # ------------------------------------------------------------------
    # Step C: CPOX[1] .. CPOX[NX-1] — general interior points
    # Each K uses a split integral: R1 (left of K) + R2 (right of K)
    # IMPORTANT: R1, R2, PS1, PS2 reset to zero at the start of each K
    # ------------------------------------------------------------------
    for K in range(1, NX):    # K corresponds to CPOX index 1..NX-1
        R1 = R2 = PS1 = PS2 = 0.0
        x0a_inner = X0 / A
        xa_inner  = x0a_inner + DELX0A

        # Left integral R1: accumulate from X0A up to point K
        for i in range(K):
            R1        += RHET_pt[i] * EVAL1(x0a_inner, xa_inner)
            x0a_inner  = xa_inner
            xa_inner  += DELX0A
        R1 *= (1.0 / x0a_inner - CTRM)

        xad_inner  = xa_inner
        xa_inner  -= DELX0A

        # Right integral R2: accumulate from point K to NX
        for i in range(K, NX):
            PS1       += RHET_pt[i + 1] * EVAL2(xa_inner, xad_inner)
            PS2       += RHET_pt[i + 1] * EVAL1(xa_inner, xad_inner)
            xa_inner   = xad_inner
            xad_inner += DELX0A
        R2 = PS1 - CTRM * PS2

        CPOX[K] = max(CI3 - A * A / DP3 * (R1 + R2), 0.0)

    # ------------------------------------------------------------------
    # Step D: CPOX[NX] — special case at X = A (right boundary / surface)
    # This uses only the right-sum integral (DO 378 in Fortran)
    # ------------------------------------------------------------------
    DXL, DXU = X0 / A, X0 / A + DELX0A
    RR2 = 0.0
    for i in range(NX):
        RR2 += RHET_pt[i] * EVAL1(DXL, DXU)
        DXL  = DXU
        DXU += DELX0A
    CPOX[NX] = max(CI3 - A * A / DP3 * (1.0 - CTRM) * RR2, 0.0)

    # ------------------------------------------------------------------
    # Step E: Update TPS from new surface concentration
    # ------------------------------------------------------------------
    DCPDX = KC3 / DP3 * (CI3 - CPOX[NX])
    H3P, DP3P = H3, DP3
    TPS   = T - (H4 * KC4 * CI4 + H3 * DP3 * DCPDX) / HC
    H3    = unbar(H3TBL, TPS)
    DP3   = DP3F(TPS, D03, P)
    TMTPO = TMTPN
    TMTPN = T - TPS

    # ------------------------------------------------------------------
    # Step F: Convergence check (skip on LP2 == 1, need two passes first)
    # ------------------------------------------------------------------
    if LP2 > 1:
        CMCPO = CMCPN
        CMCPN = CI3 - CPOX[NX]
        if (abs(TMTPO - TMTPN) / TMTPN < 0.05 and
            abs(CMCPO - CMCPN) / CMCPN < 0.05):
            break   # Phase 2 converged

    # ------------------------------------------------------------------
    # Step G: Update CPX profile for next pass (DO 55)
    # ------------------------------------------------------------------
    for i in range(NX1):
        if LP2 % 5 != 0:
            CPX[i] = 0.8 * CPX[i] + 0.2 * CPOX[i]   # normal: weighted blend
        else:
            CPX[i] = (CPOX[i] + PCPOX[i]) / 2.0       # every 5th: average with previous
        PCPOX[i] = CPOX[i]

    CMCPN = CI3 - CPX[NX]
    DCPDX = KC3 / DP3P * (CI3 - CPX[NX])
    TPS   = T - (H4 * KC4 * CI4 + H3P * DP3P * DCPDX) / HC
    TPS   = max(TPS, 1.0)
    H3    = unbar(H3TBL, TPS)
    DP3   = DP3F(TPS, D03, P)
    TMTPO = TMTPN
    TMTPN = T - TPS

    # Update reaction parameters for next pass
    GAMMA = BGM / TPS
    BETA  = -CPX[NX] * H3 * DP3 / (KP * TPS)
    K0    = ALPH2 * exp(-GAMMA) * CI1**EN3

    LP2 += 1

else:
    # LP2 hit 50 — non-convergence, return best estimate with warning
    print(f"WARNING: Phase 2 did not converge in 50 tries.")
    print(f"  CP(X/A) at surface = {CPOX[NX]:.5e}")
    print(f"  GRAD = {GRAD:.5e}   TGRAD = {TGRAD:.5e}")

# ------------------------------------------------------------------
# Final output
# ------------------------------------------------------------------
GRAD  = DCPDX * DP3
TGRAD = HC * (T - TPS)
return GRAD, TGRAD
```

---

## Known ambiguities in the scanned source — resolve these before coding

The original is a scanned document. These items are flagged as uncertain:

| Line (approx) | Issue | Most likely correct value |
|---|---|---|
| ~330 | `N = EN1` or `N = LN1` | `N = EN1` (reaction order, used as exponent in RHET) |
| ~670 | `KO = ALPH2` or `KO = ALPH3` | `ALPH2` (matches the Phase 2 usage at line ~2620) |
| ~830 | Line label `12` vs `13` | Label `12` is the TRAPP call; `13` is the DCPDX update after it |
| ~1340 | Line label `53` vs `58` | Label `53` (target of `GO TO 53` at line 1310) |
| ~1780–1860 | `CPOX(1)` vs `CPOX(I)` in DO 377 | `CPOX(1)` — this is the boundary special case, not a loop variable |
| ~1960–2010 | `DX(I)` vs `DX(1)` | `DX(I)` — inside a loop, must be indexed |

---

## Common failure modes to watch for

**Phase 1 cycling:** If TPS oscillates between two values without converging, the three-point
oscillation detector (`min/max of TPS, TPSP, TPSPP == TPSP`) should catch it and average.
If it doesn't, check that `TPSPP` and `TPSP` are being updated correctly in the right order
(TPSPP gets the old TPSP before TPSP is updated).

**LTFLG never resetting:** LTFLG is a one-shot flag. It trips to 1 when CPS < 0.25*CI3, and
resets to 0 on the very next iteration regardless of outcome. If it stays 1, you have a bug in
the reset logic.

**Negative CPS or TPS:** Both are physically impossible. The code has guards:
`CPS = max(CPS, 0.0)` and `TPS = max(TPS, 1.0)`. If CPS consistently goes negative, the
diffusivity or KC values are likely wrong (check units — the original uses °R, psia, lb, ft).

**WAF loop never breaking:** Phase 1 must `break` (converge) before execution reaches Phase 2.
If the WAF loop exhausts all retries and raises the error, the problem is upstream —
check initial guesses for CPS and DCPDX, and verify that TRAPP is returning reasonable values.

**Phase 2 CPOX profile going to zero everywhere:** This means RHET is zero, which means
either K0 or CPX are zero. On LP2==1 the CPX profile is linear from 0 at X0 to CPS at X=A —
check that CPS from Phase 1 is nonzero before entering Phase 2.

**Accumulator not resetting in Step C:** The R1, R2, PS1, PS2 accumulators in the interior
loop MUST reset to zero at the start of each K. The Fortran code does this via a GO TO that
jumps back over the initialization. If you miss this, each K accumulates all previous K values
and the profile blows up.

---

## Units

The original code uses the following unit system throughout. Do not mix with SI.

| Quantity | Unit |
|---|---|
| Temperature | °R (Rankine) |
| Pressure | psia |
| Length | ft |
| Mass flux G | lb/(ft²·s) |
| Diffusivity | ft²/s |
| Concentration | lb/ft³ or mole fraction (consistent with the COMMON block setup) |
| Reference temperature | 492 °R (= 32 °F = 273 K, the standard state) |
| Reference pressure | 14.7 psia |

The `(T/492.)^1.823` and `(14.7/P)` factors in DP3F are the Chapman-Enskog pressure and
temperature corrections from the reference state.

---

## What NOT to change

- The 5% convergence tolerance is intentional. Do not tighten it.
- The `0.8/0.2` and `0.2/0.8` weighted average coefficients are calibrated. Do not change them.
- The `NX = 24` spatial resolution in Phase 2 is deliberate.
- The every-5th-pass smoothing (`LP2 % 5`) is a deliberate stability measure. Keep it.
- The three-point oscillation detector in Phase 1 is load-bearing — do not simplify to a
  two-point check.
