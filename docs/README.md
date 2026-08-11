# Reference documents

Every empirical constant in `fvm/mechanism.py` traces to one of the reports
below. The **scanned PDFs and the multi-megabyte thermodynamic databases are
deliberately not committed** — they are ~21 MB of material that would sit in
git history forever without ever changing. What *is* committed is the Fortran,
the Python conversion work, and the small reference data files.

If you have cloned this repo and want the full set, the documents are on
[NASA NTRS](https://ntrs.nasa.gov). Drop them back into `docs/` at the paths
below and everything works, including the tests that currently skip.

## What each document supplies

### Kesten — UARL, NASA contract NAS 7-458

The hydrazine catalyst bed model comes from here.

| File (expected path) | Document |
|---|---|
| `docs/Kesten - 1967 - Analytical Study of Catalytic Reactors for Hydrazine Decomposition.pdf` | **F910461-12**, First Annual Progress Report, May 1967 |
| `docs/Kesten - 1970.pdf` | **G910461-30**, *Computer Programs Manual*, August 1968 |

> Note the second filename is misleading — the title page reads August 1968,
> not 1970. Cite it by report number.

Used by `fvm/mechanism.py`:

- Shell 405 rate constants: catalytic N₂H₄ (A = 10¹⁰ s⁻¹, Ea/R = 2500 °R),
  catalytic NH₃ (A = 0.3–1×10¹¹, Ea/R = 50,000 °R), homogeneous N₂H₄
  (A = 2.14×10¹⁰ s⁻¹, Ea/R = 33,000 °R)
- The hydrogen-inhibited ammonia rate, F910461-12 Eq. 43
- Mass-transfer correlation (subroutine `KCF`) and the diffusivity correction
- Diffusion coefficients at STP for N₂H₄ and NH₃
- Liquid hydrazine specific heat, 0.7332 Btu/lb·°R
- The dissociation-fraction convention `X = (3f + 1)/4`, F910461-12 p. 11

Two things Kesten flags about his own numbers, carried into the docstrings:
the hydrazine activation energy was *"chosen rather arbitrarily"*, and the
hydrogen inhibition order was measured on **platinum** and assumed to transfer
to Shell 405 — *"this assumption remains untested"*.

### Hayn — NASA TM-77730

| File | Document |
|---|---|
| `docs/NASA-TM-77730 - Performance Determination of Microjets.pdf` | Translation of D. Hayn, *"Beiträge zur Leistungsermittlung von Mikrodüsen"*, TU Munich dissertation, 1983. NTRS 19840026623 |

Used by `validate_hayn.py`: geometry of the 12 cold-gas micronozzles
(Table 4, p. 107), the test matrix, and the measured performance the CFD
solver is validated against.

## `kesten_claude/` — Fortran conversion work

Prior work converting Kesten's original Fortran to Python. Tracked, because it
resolves things the scans cannot:

- `fortran/` — the original listings. `PARAM.f` and `MAIN.f` contain two
  copies of the same subroutine that **disagree**: the hydrogen inhibition
  order is 1.6 in one and 1.0 in the other. That discrepancy is unresolved in
  the source material, and it changes the predicted dissociation rate by a
  factor of ~24, so it is the dominant uncertainty in the bed model.
- `kinetics.json` — parameters in machine-readable form
- `vapor_reference.csv` — **output from Kesten's original program**, used as a
  regression test in `tests/test_mechanism.py`. Our thermodynamics agrees with
  it to 1.7–1.9%.
- `kesten_verification_inputs.txt` — the full input deck for that case
- `SGRAD_codex_briefing.md` — translation notes for subroutine `SGRAD`

### A resolved question

`PARAM.f` carries the note *"This really needs to be checked. It might be
CC+HR not CC*HR"*. The **product** form in `MAIN.f` is correct: `BETA` is the
Prater number, β = (−ΔH)·D_e·C_s/(k_e·T_s), and Kesten's own Reference 8 is
Prater's paper on heat of reaction inside porous particles. The addition form
is also dimensionally impossible — you cannot add a concentration to
(enthalpy × diffusivity).

## Not committed

| Path | Size | Why it matters |
|---|---|---|
| `docs/*.pdf` | 15.5 MB | Source of every empirical constant above |
| `docs/kesten_claude/docs/Reactor_Textbook` | 3.3 MB | Background reference |
| `docs/kesten_claude/old_code/**/THERMO.INP` | 0.9 MB | NASA CEA thermo database |
| `docs/kesten_claude/old_code/**/thermo.txt` | 0.9 MB | Same, alternate format |
| `docs/kesten_claude/old_code/src/nasa9.dat` | 1.1 MB | NASA-9 polynomial data |

`THERMO.INP` is the one with a live consequence: `tests/test_chem.py`
validates the built-in NASA-7 table against it (molecular weights agree to
3 ppm, formation enthalpies to 0.042 kJ/mol). Without the file that test
**skips silently**, so a fresh clone has weaker coverage than a working copy.
Restore it if you are changing anything in `fvm/chem.py`.
