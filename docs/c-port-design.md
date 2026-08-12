# Porting the solver core to C — design note

Status: **proposal, nothing built yet.** Written to be argued with.

## 1. Scope

"Convert to C" should not mean converting the project. Measured against the
profile, almost all the runtime is in a small part of the code:

| | Lines | Goes to C? |
|---|---|---|
| `riemann.py`, `reconstruct.py`, `viscous.py` | 472 | **yes** |
| `solver.py` — residual, face fluxes, local timestep | ~150 of 350 | **yes** |
| `solver.py` — run loop, checkpointing, reporting | ~200 | no |
| grid, geometry, thermo, bc, quasi1d, post, plots | 1,600 | no |
| chemistry, mechanism, catbed, plugflow, thruster, fitting, propellants | 2,377 | no |
| tests | 1,734 | extended, not ported |

So roughly **600 lines of Python become 800–1,000 lines of C**, and 3,977 lines
stay where they are. The reacting-flow chain in particular has no business in
C: the bed integration is now 224 ms, and its cost was algorithmic, not
language.

## 2. The non-negotiable

**The Python solver stays, and the tests compare C against it.**

This is the whole risk-management strategy and everything else is detail. The
current implementation has properties that took real work to establish:

- free-stream preservation to ~10⁻¹² relative
- mass conserved throat-to-exit to 0.002%, momentum balance to 0.22%
- C_d = 0.967 against Hayn's measured micronozzle data
- the metric identity `Σ n_r S = A_planar` to round-off

If C becomes the only implementation, every one of those has to be re-earned,
and the failure mode is not a crash. A transposed index in a viscous flux
shifts thrust by a few percent and looks entirely plausible. With both
implementations present, that is a failing test instead.

Concretely: add a `backend` fixture and parameterise the existing physics tests
over `("python", "c")`. All 177 tests then protect the C path for free.

## 3. Proposed C API

Deliberately narrow — C computes, Python owns everything else.

```c
/* fvm_core.h */
typedef struct {
    int ni, nj, ng;
    const double *nx_i, *nr_i, *S_i, *L_i;    /* (ni+1, nj) */
    const double *nx_j, *nr_j, *S_j, *L_j;    /* (ni, nj+1) */
    const double *V, *A_planar, *rc;          /* (ni, nj)   */
    const double *dxi, *dri, *disti;          /* (ni-1, nj) */
    const double *dxj, *drj, *distj;          /* (ni, nj-1) */
    const double *nx_w, *nr_w, *d_wall;       /* (ni,) wall  */
    const double *Sx_I, *Sr_I, *Sx_J, *Sr_J;  /* (ni, nj) for dt */
} fvm_grid;

typedef struct { double gamma, R, Pr, mu_ref, T_mu_ref, omega; } fvm_gas;

typedef struct {
    int flux;              /* 0 roe, 1 hllc, 2 ausm            */
    int limiter;           /* 0 none, 1 minmod, ... 4 mc       */
    int viscous;           /* 0 Euler, 1 Navier-Stokes         */
    int wall_isothermal;   /* 0 adiabatic, 1 fixed T_wall      */
    double entropy_fix, T_wall;
} fvm_scheme;

/* W is ghosted primitives (4, ni+2ng, nj+2ng); R is (4, ni, nj). */
void fvm_residual(const fvm_grid *g, const fvm_gas *gas,
                  const fvm_scheme *s, const double *W, double *R);

void fvm_local_dt(const fvm_grid *g, const fvm_gas *gas, const fvm_scheme *s,
                  const double *U, double cfl, double visc_factor, double *dt);

int  fvm_num_threads(void);
void fvm_set_num_threads(int n);
const char *fvm_build_info(void);   /* compiler, OpenMP, fastmath state */
```

**Boundary conditions stay in Python.** They are where the case-specific
subtlety lives — three inlet variants, pointwise supersonic/subsonic outflow
switching, slip versus no-slip walls — and they are not hot. The cost is
crossing the boundary twice per RK stage: 5 stages × 30,000 iterations =
300,000 crossings at a few microseconds, so a second or two total against a
25-minute run. Acceptable.

If measurement later shows the crossings or the Python BC fill dominating, move
the whole step loop into C as a second phase. Do not start there.

## 4. Decisions to make deliberately

**Memory layout.** Python is variable-major, `(4, ni, nj)`. C would prefer
cell-major, `(ni, nj, 4)`, so a face reads four contiguous doubles. Options:
transpose at the API boundary (~0.1 ms per call, fine if the residual drops
from 50 ms to 10), index strided in C and let the prefetcher cope, or change
the Python layout (invasive, touches validated code — no). **Recommend:
strided first, measure, transpose only if it shows up.**

**Fast math.** `/fp:fast` or `-ffast-math` will change results at the last few
digits and can break the free-stream test, which asserts 10⁻¹². Build without
it initially so the comparison tests are meaningful; revisit as a measured
option once correctness is pinned.

**Threading.** OpenMP over the outer index of each loop nest. Face fluxes and
the cell gather are separate passes, so there are no write races if the flux
arrays are materialised. Realistic expectation, from the Numba-parallel
measurement on the same machine: **~1.4× at 220×80, ~8–9× at 880×320.** Small
grids do not parallelise — 0.6 MB across 32 threads is 20 KB each and
synchronisation dominates. First-touch allocation matters on a 32-core box.

## 5. Build on Windows

Visual Studio 2022 is installed; `cl` is not on the PATH (it needs the
Developer Command Prompt, or `vcvarsall.bat`). No CMake, gcc or clang present.

For a single DLL, CMake is not required:

```bat
call "C:\Program Files\Microsoft Visual Studio\2022\<ed>\VC\Auxiliary\Build\vcvars64.bat"
cl /LD /O2 /openmp:llvm /fp:precise csrc\fvm_core.c /Fe:fvm\_fvm_core.dll
```

Use `/openmp:llvm` rather than plain `/openmp` — the latter is OpenMP 2.0.
Adopt CMake only if this needs to build on a Linux cluster too.

**Binding: cffi in ABI mode** (`ffi.dlopen`), already installed. Cleaner
declarations than ctypes and no compile step for the binding itself. The
Python side loads the DLL if present and silently falls back to the NumPy
implementation if not, so a broken or missing build never blocks anyone.

## 6. Phasing

Each phase leaves the repo working.

1. **Comparison harness first, before any C.** Add the `backend` fixture and a
   `test_c_matches_python.py` that skips when the DLL is absent. Writing the
   test first means the port is never unverified.
2. **Inviscid path** — MUSCL, limiters, Roe. 59% of runtime, self-contained,
   and the easiest to check.
3. **Viscous path** — Green–Gauss gradients, face correction, stresses. 29%.
   The fiddly part is the wall treatment; keep the exact no-slip normal
   gradient rather than reimplementing from scratch.
4. **Local timestep** — small, but keeps `U` on the C side.
5. **OpenMP**, only once serial C matches Python. Parallelising unverified
   code is how races get blamed on physics.
6. **Optional: whole step loop in C**, only if measurement justifies it.

## 7. Honest assessment

**On speed alone, this is a weak case.** Numba with `prange` already measured
61× at 880×320 and 2.4× at 220×80, for a fraction of the effort. C with OpenMP
would likely land within 1.2–1.5× of that — both go through LLVM and do the
same fusion.

**At the current grid, neither helps much.** 220×80 is 0.6 MB, sits in cache,
and NumPy is already close to optimal there. The wins appear at grid sizes you
are not currently running.

So the strongest argument for C is not wall-clock. It is that C is reviewable,
portable, embeddable, and has a lifetime measured in decades — which matters if
this lineage is ever meant to inform flight software or be handed to people who
will not accept a Python dependency. Those are good reasons. They are just
different reasons, and worth being explicit about, because they change what
"done" looks like: an auditable, documented, standalone core rather than
whatever is fastest.

**Cheap experiment worth running first:** a Numba-parallel prototype of the
residual, half a day, gives a measured ceiling. If it lands within 20% of what
C would plausibly achieve, that number is useful either way — as an argument
for skipping the port, or as the performance target C has to beat.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Silent numerical divergence | Comparison tests from phase 1, before any C exists |
| Half-finished port, two implementations to maintain | Phase boundaries each leave the repo working; Python is never removed |
| Windows build fragility | Runtime fallback to NumPy; the DLL is an optimisation, never a requirement |
| Fast-math quietly breaking invariants | Build `/fp:precise` until correctness is pinned |
| Development slows to a crawl | Only the numerical core moves; all exploratory work stays in Python |
