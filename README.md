# Thruster Sizing

Two tools for small chemical thrusters:

1. **`thruster_sizing.py`** — algebraic sizing. Given a target thrust and nozzle
   geometry, solve for exit Mach number and chamber pressure, then report the
   usual performance parameters.
2. **`fvm/`** — an axisymmetric Navier–Stokes finite-volume solver that
   simulates the flowfield the sizing tool assumes, and reports what the
   viscous losses actually cost.

The second exists because the first cannot see the dominant loss mechanism at
this scale. At a 0.29 mm throat the throat Reynolds number is about **1200**:
the boundary layer is a large fraction of the passage, and inviscid
sizing over-predicts thrust substantially. See [Results](#results-for-the-baseline-thruster).

---

## Quick start

```
pip install -r requirements.txt

python thruster_sizing.py                 # algebraic sizing
python run_fvm_nozzle.py                  # CFD, default 200 x 80 grid
python run_fvm_nozzle.py --euler          # inviscid, for comparison
pytest -m "not slow"                       # fast verification tests
```

> `metronos` is installed directly from GitHub (see `requirements.txt`); it is
> not on PyPI. It is used only for unit handling in the driver scripts — the
> solver core is plain SI floats and NumPy.

---

## Part 1 — Algebraic sizing (`thruster_sizing.py`)

- Newton–Raphson solve for supersonic exit Mach number from the area ratio
- Newton–Raphson solve for chamber pressure from thrust, geometry and ambient
- Formatted performance summary in SI and English units
- Catalyst bed surface area and loading estimates
- Chamber wall thickness from yield strength and maximum pressure

Edit the inputs at the top of the `__main__` block and run.

### Propellant

Defaults are **LOX/LH2 at O/F = 5.0**, from thermochemical equilibrium (e.g.
NASA CEA). To change propellant, update `k`, `cstar_ideal`, `MW` and `T0`.

| Property | Value | Description |
|---|---|---|
| `k` | 1.26 | Specific heat ratio |
| `cstar_ideal` | 2350 m/s | Ideal characteristic velocity |
| `MW` | 11.8 g/mol | Combustion gas molecular weight |
| `T0` | 3250 K | Adiabatic flame temperature |
| `eta` | 0.95 | C\* efficiency |

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
| Gas | Calorically perfect, frozen composition — matches `thruster_sizing.py` |
| Viscosity | Power law `μ = μ_ref (T/T_ref)^ω`, ω = 0.7 (Sutherland also available) |
| Conductivity | From `Pr` (0.6, representative of H₂-rich products) |
| Turbulence | **None, and none is needed.** Re_throat ≈ 1200 — the flow is laminar |

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
```

---

## Verification

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

## Results for the baseline thruster

*(8.3 mm chamber, 0.29 mm throat, ε = 100, p₀ = 8.45 bar, T₀ = 3250 K,
p_amb = 5 Torr, adiabatic wall, laminar. 220 × 80 grid, Roe + van Albada,
30 000 iterations, residual down 7.2 × 10⁻⁶.)*

| Quantity | CFD | Ideal 1-D | Ratio |
|---|---|---|---|
| Thrust | 77.55 mN | 100.41 mN | **0.772** |
| Mass flow | 0.02265 g/s | 0.02434 g/s | **0.930** (C_d) |
| Specific impulse | 349.2 s | 420.7 s | **0.830** |
| Thrust coefficient | 1.389 | 1.799 | 0.772 |

**The sizing tool over-predicts thrust by 29%.** That is the headline: at this
scale the viscous loss is not a correction, it is a first-order effect, and no
choice of C\* efficiency in the algebraic tool would have revealed it.

Why:

| | |
|---|---|
| Throat Reynolds number | **1161** — laminar, thick boundary layer |
| Exit Mach, on axis | 4.37 (vs 5.41 ideal) |
| Exit Mach, first cell off wall | 0.02 |
| Boundary layer at exit plane | **1.06 mm of a 1.45 mm exit radius** |
| Sonic line at exit plane | 0.23 mm in from the wall |
| Wall axial viscous force | −42.1 mN |
| Adiabatic wall temperature | 3250 K (chamber) → 2554 K (exit) |

The boundary layer occupies about 73% of the exit radius, so most of the
nozzle's geometric area ratio is doing no useful expansion. This is the
well-known micro-nozzle result, and it is the reason a 2-D viscous solve was
worth building.

**On the C\* ratio.** The report shows C\* *above* ideal (1.074). That is not an
error: C\* is evaluated as `p₀·A_t/ṁ` with the **geometric** throat area, and the
boundary layer reduces the effective throat area, so the apparent C\* rises by
roughly `1/C_d`. The discharge coefficient is the physical quantity.

### Numerical quality of that run

| | |
|---|---|
| Mass error, inlet vs exit | 0.021% |
| Mass spread, throat to exit | 0.021% |
| Momentum balance residual | 0.020% of thrust — against a wall pressure force 588× the thrust |
| Max wall y⁺ | 0.18 |
| Grid sensitivity | 160 × 60 vs 220 × 80 agree on thrust to ~0.5% |

Mass flow and thrust are integrated from the solver's **own numerical face
fluxes**, not from cell-centred states. For a conservative scheme at steady
state, summing continuity down a column of cells makes the j-face terms
telescope to (wall flux) − (axis flux), and both vanish identically — the wall
passes no mass, and axis faces carry `S = 0`. Consecutive stations must
therefore report exactly the same mass flow. Measuring from averaged cell
states instead injects ~3% of spurious variation, enough to look like a real
leak; it also left a 1.9% momentum-balance residual that was pure measurement
error. Both drop to ~0.02% when the fluxes are integrated consistently.

Figures and a full summary land in `out/` (`nozzle_visc_RESULTS.md` plus PNGs).

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
