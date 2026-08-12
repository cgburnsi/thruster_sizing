# Propellant tables

One CSV per propellant, holding chamber conditions from a CEA rocket run.
`run_fvm_nozzle.py --propellant <name>` loads them by filename.

**Nothing is shipped for the green propellants.** Their formulations and CEA
results are yours; plausible-looking numbers invented here would be worse than
no numbers at all, because they would look authoritative. Start one with:

```python
from fvm.propellants import write_template
write_template("af_m315e")
```

then fill it from CEA. `TabulatedPropellant.from_cea("case.out")` will also
read a CEA output directly, though the CSV path is the reliable one — CEA's
format varies with the options used.

## Format

| column | required | notes |
|---|---|---|
| `p_bar` | yes | chamber pressure |
| `T_K` | yes | chamber temperature |
| `MW` | yes | kg/kmol |
| `gamma` | yes | ratio of specific heats |
| `cstar_m_s` | no | **prefer CEA's value** — it accounts for shifting equilibrium through the throat, which the closed form does not |
| `mu_Pa_s`, `Pr` | no | transport; omitted means a power law anchored at `T_K` |
| `x_<species>` | no | mole fractions, carried for reference |

One row makes the propellant pressure-independent, which is a reasonable
approximation over a modest range. Several rows are interpolated log-linearly
in pressure; going outside the tabulated range raises rather than
extrapolating.

## `hydrazine_X084.csv`

Generated from `fvm.mechanism.HydrazineShell405` at X = 0.84, **not** from CEA.
It exists so the tabulated path can be checked against the reacting model
through the same interface, and as a worked example of the format.
