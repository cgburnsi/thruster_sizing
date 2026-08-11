# Thruster Sizing

Modelling tools for a small **catalytic monopropellant thruster** — liquid
hydrazine over a Shell 405 bed, expanded through a 0.29 mm throat.

1. **`thruster_sizing.py`** — algebraic sizing. Given a target thrust and nozzle
   geometry, solve for exit Mach number and chamber pressure, then report the
   usual performance parameters.
2. **`fvm/`** — an axisymmetric Navier–Stokes finite-volume solver that
   simulates the flowfield the sizing tool assumes, and reports what the
   viscous losses actually cost.
3. **`fvm/` (reacting)** — a catalyst bed model: multi-species thermodynamics,
   a propellant-agnostic reaction mechanism, bed packing and pressure drop, and
   a 1-D two-temperature reacting plug-flow reactor, coupled to the nozzle so
   the chain runs from liquid feed to thrust.

The CFD solver exists because algebraic sizing cannot see the dominant loss
mechanism at this scale. At a 0.29 mm throat the boundary layer occupies a
large fraction of the passage and inviscid sizing over-predicts thrust. How
badly depends strongly on the propellant, because that sets the Reynolds
number. At **Re = 6,184** for catalytic hydrazine the over-prediction is 9%;
the same solver gives 29% for LOX/LH2 at Re = 1,164. See
[Results](#results-for-the-baseline-thruster).

The solver is [verified](#verification) against analytic solutions and
[validated](#validation-against-measurement) against measured micronozzle data
from NASA TM-77730, whose test hardware matches this thruster's geometry to
within 2%.

---

## Quick start

```
pip install -r requirements.txt

python thruster_sizing.py                 # algebraic sizing
python run_fvm_nozzle.py                  # CFD, default 200 x 80 grid
python run_fvm_nozzle.py --euler          # inviscid, for comparison
python run_thruster.py                    # bed + nozzle coupled, liquid feed to thrust
python validate_hayn.py --pc 10           # validation vs NASA TM-77730
pytest -m "not slow"                       # fast verification tests
```

> `metronos` is installed directly from GitHub (see `requirements.txt`); it is
> not on PyPI. It is used only for unit handling in the driver scripts — the
> solver core is plain SI floats and NumPy.

A production CFD run takes tens of minutes on one core. Both drivers checkpoint
periodically and resume automatically, so an interrupted run costs at most
`--checkpoint-every` iterations; `--fresh` forces a clean start.

---

## Part 1 — Algebraic sizing (`thruster_sizing.py`)

- Newton–Raphson solve for supersonic exit Mach number from the area ratio
- Newton–Raphson solve for chamber pressure from thrust, geometry and ambient
- Formatted performance summary in SI and English units
- Catalyst bed surface area and loading estimates
- Chamber wall thickness from yield strength and maximum pressure

Edit the inputs at the top of the `__main__` block and run.

### Propellant

Gas properties are **not hardcoded**. They are pulled from
`fvm.mechanism.HydrazineShell405` at a chosen ammonia dissociation fraction, so
this tool and the bed model cannot drift apart:

```python
X_dissociation = 0.84                     # what fvm.thruster predicts for this bed
_chamber = HydrazineShell405().chamber_conditions(X_dissociation)
k, MW, T0, cstar_ideal = ...              # all follow
```

X is the single design variable — it sets T₀, molecular weight and γ together.
Run `run_thruster.py` to recompute it for a different bed or flow rather than
guessing. At X = 0.84:

| Property | Value |
|---|---|
| `k` | 1.337 |
| `cstar_ideal` | 1258 m/s |
| `MW` | 11.49 g/mol |
| `T0` | 995 K |
| `eta` | 0.95 (empirical C\* efficiency) |

The catalyst bed section uses `fvm.catbed`, and reproduces the original hand
calculation exactly (a_v = 5423 1/m, 58.68 cm² for 60 PPI foam).

---

## Part 2 — Axisymmetric Navier–Stokes solver (`fvm/`)

A second-order cell-centred finite-volume solver on a structured, body-fitted
grid. Vectorised NumPy throughout; no compiled extensions.

### Governing equations

Cylindrical `(x, r)` with azimuthal symmetry. Multiplying through by `r` puts
the equations in strong conservation form on the true axisymmetric control
volume:

```
d(rU)/dt + d(r·Fx)/dx + d(r·Fr)/dr = [0, 0, p − τ_θθ, 0]
```

with `U = [ρ, ρu, ρv, ρE]`. Discretely, over a cell:

```
V · dU/dt = − Σ_faces (F_inv − F_visc)·n S  +  [0, 0, p − τ_θθ, 0] · A_planar
```

where `V = ∫ r dA`, `S = ∫ r dl`, and `A_planar` is the plain cell area.

Two consequences of this form are worth stating, because both are easy to get
wrong and both are covered by tests:

- **The pressure source is integrated over the planar area, not the r-weighted
  volume.** It is consistent only because the metrics satisfy
  `Σ n_r S = A_planar` exactly, which is asserted to round-off in
  `test_closed_surface_identities`.
- **The symmetry axis needs no special flux treatment.** Faces on `r = 0` have
  `S = 0`, so their flux contribution vanishes identically. Ghost cells at the
  axis exist only to supply gradients.

### Physical model

| | |
|---|---|
| Gas | Calorically perfect, frozen composition — supplied by the bed model, or set directly |
| Viscosity | Power law `μ = μ_ref (T/T_ref)^ω`, ω = 0.7 (Sutherland also available) |
| Conductivity | From `Pr` (0.6–0.7, representative of H₂-rich products) |
| Turbulence | **None, and none is needed.** Re_throat is 1,000–6,000 — the flow is laminar |

### Numerics

| Component | Options | Default |
|---|---|---|
| Inviscid flux | Roe (Harten–Hyman entropy fix), HLLC, AUSM⁺-up | `roe` |
| Reconstruction | MUSCL on primitives, with positivity fallback | 2nd order |
| Limiter | van Albada, minmod, van Leer, MC, none | `vanalbada` |
| Viscous flux | Green–Gauss cell gradients + directional face correction | — |
| Time integration | 1/3/5-stage low-storage Runge–Kutta | `rk5` |
| Convergence | Local time stepping, CFL ramp, optional residual smoothing | — |

The entropy fix is applied to the **acoustic waves only** — smearing the
entropy and shear waves would thicken the boundary layer, which is the thing
being measured.

### Boundary conditions

| Boundary | Treatment |
|---|---|
| Inlet | Subsonic, imposes `p0`/`T0`, extrapolates area-averaged velocity |
| Outlet | Supersonic extrapolation, switching **pointwise** to characteristic back-pressure wherever local `M < 1` (which it is, inside the exit-plane boundary layer) |
| Axis | Symmetry (mirror `v`) |
| Wall | No slip (viscous) or tangency (Euler); adiabatic or isothermal |

**On the inlet.** The contraction ratio here is 819:1, which puts the chamber at
M ≈ 7×10⁻⁴ and a dynamic pressure of ~0.3 Pa against an 845 kPa static level.
Both textbook subsonic inlets fail in that regime, and both are still available
so the failure is reproducible:

- `riemann` extrapolates the outgoing `qn − a` invariant. It recovers a ~1 m/s
  chamber velocity by differencing two quantities of order 13 000 m/s — no
  significant digits survive.
- `pressure` extrapolates static pressure and reads M off `p0/p`. Since
  `u ∼ √(p0 − p)`, sub-Pascal cell-to-cell noise becomes order-100% scatter in
  inlet velocity, which seeds a spurious axial jet.
- `stagnation` (default) extrapolates the **area-averaged normal velocity** and
  applies it uniformly. Velocity is an O(1 m/s) quantity carried across at face
  value, so nothing is amplified. Switching to it took the residual stall from
  4×10⁻² down to 6×10⁻⁴ and the station-to-station mass-flow spread from 254%
  to 8.7%.

### Module map

| Module | Contents |
|---|---|
| `fvm/thermo.py` | `PerfectGas`: equation of state, transport, isentropic relations |
| `fvm/geometry.py` | `ConicalNozzle`, `BellNozzle`, `NozzleContour.from_points` |
| `fvm/grid.py` | Algebraic grid generation, clustering, all FV metrics |
| `fvm/riemann.py` | Roe / HLLC / AUSM⁺-up flux functions |
| `fvm/reconstruct.py` | Limiters and MUSCL reconstruction |
| `fvm/viscous.py` | Gradients, viscous fluxes, hoop-stress source |
| `fvm/bc.py` | Ghost-cell boundary conditions |
| `fvm/solver.py` | `NozzleSolver`: residual assembly, time stepping, checkpointing |
| `fvm/quasi1d.py` | Isentropic reference solution and field initialisation |
| `fvm/post.py` | Performance integrals, wall quantities, formatted report |
| `fvm/plots.py` | Figures |

### Library use

```python
from fvm import (PerfectGas, LOX_LH2_OF5, ConicalNozzle, Grid,
                 BoundaryConditions, NozzleSolver, quasi1d, post)

gas     = PerfectGas(**LOX_LH2_OF5)
contour = ConicalNozzle(r_chamber=4.15e-3, r_throat=0.145e-3,
                        area_ratio=100.0, L_chamber=4e-3)
grid    = Grid.from_contour(contour, ni=220, nj=80, wall_spacing=0.003)
bcs     = BoundaryConditions(p0=845e3, T0=3250.0, p_amb=666.6)

sol = NozzleSolver(grid, gas, bcs)
sol.initialize(quasi1d.initial_field(grid, gas, bcs.p0, bcs.T0, bcs.p_amb))
sol.run(max_iter=30000, tol=1e-9)

post.report(sol, quasi1d.ideal_performance(contour, gas, bcs.p0, bcs.T0, bcs.p_amb))
```

### Driver options

```
python run_fvm_nozzle.py [options]

  --ni / --nj N          grid size (default 200 x 80)
  --wall-spacing F       first cell height off the wall, as a fraction of local radius
  --l-chamber MM         override the modelled chamber length
  --iters N --tol T      iteration budget and convergence target
  --flux {roe,hllc,ausm} --limiter {...} --integrator {rk1,rk3,rk5}
  --euler                inviscid run
  --bell                 bell contour instead of conical
  --wall-temp K          isothermal wall (default adiabatic)
  --restart FILE.npz     resume from a checkpoint
  --no-plots

  --chamber P0_BAR T0_K MW GAMMA
                         chamber conditions, as printed by run_thruster.py.
                         Without it the built-in LOX/LH2 block is used, which
                         is not what this thruster burns.
  --pr --mu-ref --omega  transport properties for the --chamber gas
```

### Closing the loop

`run_thruster.py` prints a ready-to-paste command for the CFD run:

```
python run_fvm_nozzle.py --chamber 8.4400 994.2 11.492 1.3374     --ni 220 --nj 80 --wall-spacing 0.003 --l-chamber 4.0 --iters 30000
```

Take the discharge coefficient that run reports and hand it back:

```
python run_thruster.py --mdot 0.0443 --cd 0.93
```

which re-converges the bed and nozzle against the real throat. As a check, the
quasi-1-D reference `run_fvm_nozzle.py` prints under `--chamber` reproduces the
coupled system's own answer independently — 0.04430 g/s, 95.97 mN, Isp 220.9 s
from both.

---

## Part 3 — Catalyst bed and system coupling

The nozzle solver takes chamber conditions as given. This part produces them,
from liquid hydrazine.

### Module map

| Module | Contents |
|---|---|
| `fvm/chem.py` | NASA-7 species thermodynamics, mixtures, `h(T)` inversion |
| `fvm/mechanism.py` | Propellant-agnostic `Mechanism`; `HydrazineShell405` |
| `fvm/catbed.py` | Bed geometry, packing, Ergun pressure drop, surface area |
| `fvm/plugflow.py` | 1-D two-temperature reacting plug-flow reactor |
| `fvm/thruster.py` | Couples bed to nozzle; solves the operating point |

### The mechanism

Three paths, each per mole of its limiting reactant:

```
N2H4 -> 4/3 NH3 + 1/3 N2      catalytic, diffusion-controlled
N2H4 -> 4/3 NH3 + 1/3 N2      homogeneous, gas phase
NH3  -> 1/2 N2  + 3/2 H2      catalytic, kinetically controlled,
                              inhibited by hydrogen
```

The competition between the first and last is the whole design problem.
Decomposition sets the temperature; dissociation eats some of it back while
lowering molecular weight, so **c\* peaks at partial dissociation** — the model
puts the optimum at X ≈ 0.28, and hydrazine engines are designed for 30–50%.

Two structural choices:

- **Heats of reaction are computed from species enthalpies, never tabulated.**
  A hard-coded ΔH can silently disagree with the thermodynamic data beside it.
- **Rates combine kinetics and mass transfer as resistances in series.** This
  is physically right for a catalytic bed and it puts the answer where the data
  is trustworthy: hydrazine decomposition over iridium is diffusion-limited, so
  the Sherwood correlation controls rather than Arrhenius constants. A test
  asserts it — raising the pre-exponential 1000× must not move the rate.

### The reactor

```
G dY_i/dz = omega_i        species
G dh/dz   = -q_loss        energy   (adiabatic => h constant)
  dp/dz   = -(Ergun)       momentum
```

Gas and catalyst carry separate temperatures. Reactions see the solid; the gas
is heated by convection from it. The solid temperature is not integrated but
solved from a local balance at each station, and it is genuinely two-sided —
the catalyst runs **750 K above the gas** where decomposition dominates and
below it where dissociation does.

Note the energy equation has no reaction source term: for an adiabatic bed the
reaction enthalpy is already in the species enthalpies, so total enthalpy is
conserved. That gives a strong check, and a test takes it — the integrated exit
temperature matches the closed-form adiabatic value to 1e-6.

**Vaporisation is not modelled.** Integration begins at the vapour-region
inlet. That is a necessity, not a shortcut: conserving enthalpy through
instantaneous vaporisation puts the vapour *below* the feed temperature,
Arrhenius rates vanish and nothing ignites. Real beds resolve this by
conducting heat upstream, which a steady 1-D model cannot represent.

The inlet composition is not a free parameter either. The bed must decompose
some of its own feed to vaporise the rest, and that extent follows from
requiring `h(T_vapor, Y) = h_liquid(T_feed)`. It comes out at **34.2%
decomposed at 455.6 K, against 38% in Kesten's own vapour region** — nothing
in that calculation is fitted to his data.

### Coupling

Bed and nozzle need each other: chamber pressure is the bed exit, the choked
throat sets mass flow from chamber pressure and c\*, and the bed's pressure
drop depends on the resulting mass flux. Solved as a fixed point in both
directions — feed pressure for a given flow, or flow for a given feed pressure.

The loop uses the quasi-1-D nozzle because it needs many bed integrations. Run
the CFD **once** afterwards on the converged chamber conditions, then feed its
discharge coefficient back through `--cd`.

### Provenance

Rate parameters and correlations come from Kesten's UARL work under NASA
contract NAS 7-458 (see `docs/README.md`). He is candid about what they are,
and so is the code: the hydrazine activation energy was *"chosen rather
arbitrarily"*, and the hydrogen inhibition order was measured on **platinum**
and assumed to transfer to Shell 405 — *"this assumption remains untested"*.

**Treat them as defaults to re-fit against your own bed data, not as physical
constants.** In order of how much they move the answer:

| Parameter | Uncertainty | Effect |
|---|---|---|
| `A_NH3` | factor of 3 (Kesten's own range) | sets dissociation |
| `n_H2` | 1.0 vs 1.6, unresolved in his own Fortran | **factor of 24** in rate |
| `T_vapor` | modelling input | sets pre-decomposition |
| bed voidage | assumed | 20–35% systematic in Δp |

### Validation

Kesten's own program output (`docs/kesten_claude/vapor_reference.csv`) is a
regression test:

| | Kesten | Ours |
|---|---|---|
| Exit temperature | 1058.5 K | within 2% |
| Dissociation X | 0.7341 | within 0.03 |
| Adiabatic curve, X = 0 / X = 1 | ~1650 / ~880 K | 1646 / 868 K |
| Dissociation fraction, all 6 stations | — | matches to 5 decimals |

Sweeping the two uncertain parameters, the combination reproducing his output
is `n_H2 = 1.0, A_NH3 = 1e11` — exactly `MAIN.f` plus his input deck, not
`PARAM.f`'s 1.6. That is evidence about which of his two conflicting code
paths ran, not proof that 1.0 is better physics: this model omits the
intraparticle diffusion resistance his carries.

---

## Verification

*Verification asks whether the equations are being solved correctly.
[Validation](#validation-against-measurement) — a separate section below —
asks whether they are the right equations.*

`pytest` — 40 tests, layered so a failure points at a layer rather than at
"the CFD is wrong". Whole-solver runs are marked `slow`; use
`pytest -m "not slow"` for the fast set.

The two sharpest checks:

- **Free-stream preservation.** A uniform axial state must produce exactly zero
  residual. It holds to ~10⁻¹² relative, and only holds if the metric
  identities close *and* the axisymmetric pressure source is integrated over
  the planar area.
- **Inviscid limit.** With viscosity off, the solver recovers the quasi-1-D
  isentropic mass flow to within 2% and exit Mach number to within 2% on a
  moderate grid, with mass conserved to <1% through the divergent section.

Also covered: flux consistency `F(W,W,n) = F_phys`, conservation symmetry
`F(L,R,n) = −F(R,L,−n)`, full-upwind behaviour in supersonic flow, limiter
TVD properties, Green–Gauss exactness on linear fields, grid volume against
analytic integrals, and the slip/no-slip wall distinction.

---

## Validation against measurement

`validate_hayn.py` runs the solver against Dieter Hayn's micronozzle
experiments — *"Beiträge zur Leistungsermittlung von Mikrodüsen"*, TU Munich
dissertation 1983, translated as **NASA TM-77730** (in `docs/`). Twelve
cold-gas N₂ micronozzles, each machined down through five area-ratio steps and
fired at five chamber pressures: 1200 individual tests.

**His nozzle 2 at its ε = 100 step is this thruster to within ~2%:**

| | Hayn #2 | This thruster |
|---|---|---|
| Throat radius | 0.142 mm | 0.145 mm |
| Exit diameter | 2.84 mm | 2.90 mm |
| Area ratio | 99.65 | 100.0 |
| Divergence half-angle | 15.0° | 15° |

Geometry comes from his Table 4 (p. 107), which gives throat radius, exit
diameter, and an overall length from an unstated datum. Differencing length
against diameter across all five ε steps recovers the 15.0° divergent angle and
locates the throat 8.45 mm from that datum; subtracting the 5.0 mm chamber
fixes the convergent half-angle at ~28°. The reconstruction reproduces the
tabulated lengths at ε = 200 and ε = 10 to within 0.03 mm.

### Results

Cold N₂, T₀ = 293 K, p_amb = 1 Torr, 200 × 100 grid, 22 000 iterations:

| | 2 bar | 10 bar |
|---|---|---|
| Re_throat | 8,513 | 42,568 |
| **Discharge coefficient C_d** | **0.972** | **0.987** |
| **Specific impulse** | **70.8 s** | **74.4 s** |
| Isp / ideal | 0.958 | 0.975 |
| Thrust / ideal | 0.932 | 0.963 |
| Boundary layer at exit | 0.545 mm | 0.281 mm |
| Max wall y⁺ | 0.38 | 0.67 |
| Mass spread, throat→exit | 0.027% | 0.009% |
| Momentum balance residual | 1.16% | 0.35% |

Both fall inside Hayn's measured **63–75 s** band for cold-N₂ micronozzles at
this scale (Fig. 73), and C_d is below unity in both — a viscous nozzle cannot
pass more mass than ideal, so that is a floor the solver had to clear on its
own.

**The trend is the stronger evidence.** Dropping Reynolds number 5× makes every
loss grow, in the right direction and by a plausible amount: C_d 0.987 → 0.972,
Isp ratio 0.975 → 0.958, exit boundary layer 0.281 → 0.545 mm. Hitting one
number inside a band is easy; reproducing the Reynolds scaling is the part that
suggests the physics is being solved rather than fitted.

### What this does and does not establish

- **It does not validate the thruster prediction.** Cold N₂ runs Re = 8,500 to
  43,000; the LOX/LH2 thruster runs **1,160** — 7× below the lowest point tested
  here. Same geometry, different viscous regime. The trend extrapolates the
  right way, but it is extrapolation.
- **A laminar solver is nonetheless correct here.** Hayn addresses this directly
  (p. 31): the favourable pressure gradient puts the Boldman acceleration
  parameter several times above critical, relaminarising the boundary layer, so
  his own boundary-layer model assumes a laminar profile too.
- **The comparison band is soft.** 63–75 s is read off a 1983 scanned pen plot,
  and that figure is labelled for a different nozzle/ε combination than the one
  modelled. It brackets the right regime; it is not a point-to-point match.
- **Neither run is fully converged.** Residual drop 1.4×10⁻⁴ and 1.7×10⁻⁴,
  against 7×10⁻⁶ for the hot case, with ~1% inlet-to-exit mass error remaining.
  The 10 bar residual oscillated mid-run before coming down. Treat the third
  significant figure as unsettled.

### Independent corroboration of the mechanism

Hayn found all 12 nozzles at ε > 100 produced **more** thrust than his model
predicted, because the model had assumed Summerfield separation (separation once
wall pressure falls below 0.4 × ambient). His explanation (p. 117): *"due to the
greater influence of boundary layer formation in the case of extremely small
propulsion unit dimensions... the residual pressure near the wall is higher than
in purely potential flow."*

The solver reproduces exactly that on the hot thruster, without being told to:

| | |
|---|---|
| Inviscid 1-D exit pressure | 419 Pa → p_e/p_a = 0.63 (**over**-expanded) |
| CFD wall exit pressure | 831 Pa (2× the inviscid value) |
| CFD area-averaged exit pressure | 1419 Pa → p_e/p_a = 2.13 (**under**-expanded) |
| Wall pressure below 0.4·p_amb | never |
| Reverse flow / separation | none; wall shear positive throughout |

Boundary-layer blockage flips the nozzle from over-expanded to under-expanded —
the same mechanism Hayn measured, arrived at independently.

He also reports optimized bell contours giving **up to 2% higher Isp** than
conical at this scale; `BellNozzle` is implemented if you want to test that.

### Next step

A `--pc 5` run fills in the middle of the Reynolds range and shows whether the
trend is smooth rather than two points and a hopeful line:

```
python validate_hayn.py --pc 5 --p-amb 1.0 --ni 200 --nj 100 \
    --wall-spacing 0.0011 --iters 30000
```

---

## Results for the baseline thruster

*Catalytic hydrazine, 8.3 mm chamber, 0.29 mm throat, ε = 100, p_amb = 5 Torr,
adiabatic wall, laminar. Chamber conditions from `run_thruster.py`; nozzle from
`run_fvm_nozzle.py` on a 220 × 80 grid, 30 000 iterations, residual down
1.3 × 10⁻⁵.*

| Quantity | CFD | Ideal 1-D | Ratio |
|---|---|---|---|
| Thrust | **87.42 mN** | 95.97 mN | **0.911** |
| Mass flow | 0.04285 g/s | 0.04430 g/s | **0.967** (C_d) |
| Specific impulse | **208.1 s** | 220.9 s | **0.942** |
| Thrust coefficient | 1.568 | 1.722 | 0.911 |

**Inviscid sizing over-predicts thrust by 9%.** Real, but far less than the 29%
the same solver gives for LOX/LH2 — because propellant sets Reynolds number:

| | Re_throat | Boundary layer at exit | Thrust ratio |
|---|---|---|---|
| LOX/LH2, T₀ 3250 K | 1,164 | 1.06 mm of 1.45 mm (73%) | 0.772 |
| **Hydrazine, T₀ 994 K** | **6,184** | **0.60 mm of 1.45 mm (42%)** | **0.911** |

Cooler gas is denser and less viscous, so the layer thins and most of the
geometric area ratio does useful work. Exit Mach still runs 5.43 on the axis
against 0.015 at the wall.

The C\* ratio comes out **above** unity at 1.034. That is not an error: C\* uses
the geometric throat area while the boundary layer shrinks the effective one,
so it rises as 1/C_d. Check: 1/0.9672 = 1.0339, matching to four figures.

### The loop closes

Feeding C_d = 0.9672 back through `run_thruster.py --cd` raises chamber
pressure from 8.440 to 8.726 bar — exactly the 1/C_d = 1.0339 the throat
demands, with dissociation and bed temperature unmoved.

### Design consequence

The 0.0443 g/s point delivers **87.4 mN, not 100**. Thrust scales linearly with
chamber pressure at fixed geometry, so:

| | Current | For 100 mN |
|---|---|---|
| Mass flow | 0.0443 g/s | **0.0507 g/s** |
| Feed pressure | 8.46 bar | **10.01 bar** |
| Chamber pressure | 8.44 bar | 9.99 bar |
| Real thrust | 87.4 mN | 100.8 mN |
| Real Isp | 208.1 s | 202.9 s |

Isp falls slightly because the higher chamber pressure shifts dissociation from
X = 0.841 to 0.836.

### Numerical quality

| | |
|---|---|
| Mass spread, throat to exit | 0.0024% |
| Mass error, inlet vs exit | 0.031% |
| Momentum balance residual | 0.221% |
| Max wall y⁺ | 0.61 |

Mass flow and thrust are integrated from the solver's **own numerical face
fluxes**, not from cell-centred states. For a conservative scheme at steady
state the j-face terms telescope to (wall flux) − (axis flux) and both vanish,
so consecutive stations must report exactly the same mass flow. Measuring from
averaged cell states instead injects ~3% of spurious variation and leaves a
1.9% momentum-balance residual that is pure measurement error.

---

## Known limitations

- **Frozen, calorically perfect chemistry.** Real LOX/LH2 recombines through
  the nozzle, which recovers some Isp. Thermally perfect (`cp(T)`) and
  finite-rate models are the natural next steps.
- **Continuum assumption.** At the exit plane of an ε = 100 micro-nozzle the
  Knudsen number is no longer negligible, so the far downstream end of the
  expansion is at the edge of where Navier–Stokes applies. The throat and the
  early divergent section — where thrust is set — are safely continuum.
- **The chamber converges slowly.** At CR = 819 it sits near 1 m/s, so its
  convective time scale is roughly 800× the nozzle's. Nozzle-region mass
  conservation settles long before the global inlet-to-exit balance does; use
  `--l-chamber` to shrink the reservoir if the global figure matters to you.
- **Total enthalpy is not discretely conserved.** The scheme conserves mass,
  momentum and total energy; total enthalpy drifts up by ~4.6% at worst near the
  axis at the exit plane, growing through the divergent section. This is normal
  for a second-order upwind scheme in a strong expansion on a stretched grid,
  and it shrinks with grid refinement — but it means exit velocity carries
  roughly a 2% scheme-level uncertainty on top of everything else.
- **Explicit time stepping.** Adequate here (~90 min on one core for the
  30 000-iteration production run), but an
  implicit scheme would be needed for much finer grids or for cooled-wall cases
  with very small near-wall cells.
- **AUSM⁺-up is not recommended for this geometry.** It is unstable in the
  near-stagnant chamber for `M_inf` below ~0.2, and dissipates enough above
  that to cost ~2% in predicted mass flow relative to Roe.
- **No measured data exists at the thruster's Reynolds number.** The
  [validation](#validation-against-measurement) covers Re = 8,500–43,000 in cold
  N₂; the LOX/LH2 thruster runs at Re ≈ 1,160. Nothing here confirms the
  solver at that condition — the 77.6 mN prediction rests on verified numerics
  and a correctly-signed Reynolds trend, not on measurement.
