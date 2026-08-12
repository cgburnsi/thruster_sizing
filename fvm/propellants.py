"""Propellants defined by tabulated chamber conditions.

Why tabulated rather than modelled
----------------------------------
For hydrazine this code derives chamber conditions from a reaction mechanism,
because the whole design question is *how far* ammonia dissociates and that is
one variable. Green propellants do not work that way: HAN and ADN products sit
near equilibrium across many species, so there is no single extent to solve
for, and the composition follows from an equilibrium calculation.

The right tool for that is NASA CEA, which you already have and which is
authoritative. So rather than reimplement it badly, this module imports its
answer. A propellant becomes a small table of ``(p, T0, MW, gamma, c*)``,
presented through the same interface :class:`~fvm.mechanism.HydrazineShell405`
offers, so the nozzle solver cannot tell the difference.

That has a useful consequence: **any** propellant works, including ones this
code knows nothing about. Nothing here estimates a property.

Nothing is shipped for the green propellants themselves. Their formulations
and CEA results are yours to supply -- inventing plausible-looking numbers for
AF-M315E or LMP-103S is precisely the error this module exists to avoid. Run
``write_template()`` for a starting file, or point :meth:`from_cea` at a CEA
output.

Pressure dependence
-------------------
Give one row and the propellant is treated as pressure-independent, which is
a decent approximation over a modest range: T0 and MW move only weakly with
chamber pressure for these propellants. Give several and properties are
interpolated, log-linearly in pressure. Extrapolation past the tabulated range
raises rather than guessing.
"""
import csv
import re
from pathlib import Path

import numpy as np

from .thermo import R_UNIVERSAL as RU, PerfectGas

#: Where built-in and user propellant tables live.
PROPELLANT_DIR = Path(__file__).resolve().parents[1] / "propellants"

REQUIRED = ("p_bar", "T_K", "MW", "gamma")


class TabulatedPropellant:
    """Chamber conditions from a table, usually a CEA run.

    Parameters
    ----------
    name : str
    p, T, MW, gamma : array_like
        Chamber pressure [Pa], temperature [K], molecular weight [kg/kmol] and
        ratio of specific heats, one entry per tabulated pressure.
    cstar : array_like, optional
        Characteristic velocity [m/s]. Computed from T, MW and gamma when
        absent -- but prefer CEA's value if you have it, since CEA accounts for
        shifting equilibrium through the throat and the closed form does not.
    mu, Pr : array_like or float, optional
        Transport properties. Defaults follow a power law anchored at T0.
    composition : dict, optional
        Mole fractions, carried for reference and reporting only.
    """

    def __init__(self, name, p, T, MW, gamma, cstar=None, mu=None, Pr=0.7,
                 composition=None, source=""):
        self.name = name
        self.p = np.atleast_1d(np.asarray(p, dtype=float))
        order = np.argsort(self.p)
        self.p = self.p[order]

        def col(v):
            if v is None:
                return None
            a = np.atleast_1d(np.asarray(v, dtype=float))
            if a.size == 1:
                a = np.full(self.p.shape, a[0])
            if a.shape != self.p.shape:
                raise ValueError(f"{name}: expected {self.p.size} values, got {a.size}")
            return a[order]

        self.T = col(T)
        self.MW = col(MW)
        self.gamma = col(gamma)
        self.cstar_table = col(cstar)
        self.mu_table = col(mu)
        self.Pr = float(Pr)
        self.composition = dict(composition or {})
        self.source = source

    # -- interpolation ----------------------------------------------------
    def _interp(self, arr, p0):
        if arr is None:
            return None
        if self.p.size == 1:
            return float(arr[0])
        lo, hi = self.p[0], self.p[-1]
        if not (lo * (1 - 1e-9) <= p0 <= hi * (1 + 1e-9)):
            raise ValueError(
                f"{self.name}: chamber pressure {p0 / 1e5:.3f} bar is outside "
                f"the tabulated range {lo / 1e5:.3f}-{hi / 1e5:.3f} bar. "
                f"Extend the table rather than extrapolating.")
        return float(np.interp(np.log(p0), np.log(self.p), arr))

    @staticmethod
    def _cstar(gamma, R, T):
        g = gamma
        return (np.sqrt(g * R * T)
                / (g * np.sqrt((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0)))))

    def chamber_conditions(self, p0=None):
        """(T, MW, gamma, R, c*) at a chamber pressure, matching the hydrazine
        mechanism's interface so the nozzle side is propellant-agnostic."""
        p0 = float(self.p[0] if p0 is None else p0)
        T = self._interp(self.T, p0)
        MW = self._interp(self.MW, p0)
        g = self._interp(self.gamma, p0)
        R = RU / MW
        cs = self._interp(self.cstar_table, p0)
        return dict(T=T, p=p0, MW=MW, gamma=g, R=R,
                    cstar=cs if cs is not None else float(self._cstar(g, R, T)),
                    composition=self.composition, name=self.name)

    def perfect_gas(self, p0=None, omega=0.7):
        """A frozen-composition :class:`PerfectGas` for the nozzle solver."""
        c = self.chamber_conditions(p0)
        mu = self._interp(self.mu_table, c["p"])
        if mu is None:
            # Same anchor the reacting model uses, so a tabulated propellant
            # and a modelled one describe transport the same way.
            mu = 3.5e-5 * (c["T"] / 1000.0) ** omega
        return PerfectGas(gamma=c["gamma"], MW=c["MW"], Pr=self.Pr,
                          mu_ref=mu, T_mu_ref=c["T"], mu_law="power", omega=omega)

    # -- loading ----------------------------------------------------------
    @classmethod
    def from_csv(cls, path, name=None):
        path = Path(path)
        with path.open(newline="", encoding="utf-8-sig") as fh:
            rows = [r for r in csv.DictReader(
                ln for ln in fh if not ln.lstrip().startswith("#"))]
        if not rows:
            raise ValueError(f"{path}: no data rows")
        header = {k.strip().lower(): k for k in rows[0]}
        # Compare case-insensitively: the canonical names are mixed case for
        # readability but files in the wild will not be.
        missing = [c for c in REQUIRED if c.lower() not in header]
        if missing:
            raise ValueError(
                f"{path}: missing required column(s) {missing}. "
                f"Required: {list(REQUIRED)}. Found: {list(rows[0])}")

        def grab(col, required=True):
            key = header.get(col.lower())
            if key is None:
                return None
            vals = [r[key].strip() for r in rows]
            if any(v == "" for v in vals):
                if required:
                    raise ValueError(f"{path}: blank value in column {col!r}")
                return None
            return [float(v) for v in vals]

        comp = {k[2:]: float(rows[0][v]) for k, v in header.items()
                if k.startswith("x_") and rows[0][v].strip()}
        return cls(name or path.stem,
                   p=[v * 1e5 for v in grab("p_bar")],
                   T=grab("T_K"), MW=grab("MW"), gamma=grab("gamma"),
                   cstar=grab("cstar_m_s", required=False),
                   mu=grab("mu_pa_s", required=False),
                   Pr=float(rows[0][header["pr"]]) if "pr" in header
                   and rows[0][header["pr"]].strip() else 0.7,
                   composition=comp, source=str(path))

    @classmethod
    def from_cea(cls, path, name=None):
        """Best-effort parse of a CEA rocket-problem output file.

        CEA's output format shifts with the options used, so this reads the
        chamber column of the common ``rocket`` case and raises with what it
        did find if that fails. When it does, fall back to a CSV -- four
        numbers typed by hand are more reliable than a fragile parser.
        """
        path = Path(path)
        text = path.read_text(encoding="utf-8", errors="replace")

        def first_column(label, *, after=None):
            body = text.split(after, 1)[1] if after and after in text else text
            m = re.search(rf"^\s*{label}\s+([-\d.Ee+]+)", body, re.M)
            return float(m.group(1)) if m else None

        p_bar = first_column(r"P,\s*BAR")
        T = first_column(r"T,\s*K")
        MW = first_column(r"M,\s*\(1/n\)") or first_column(r"MW,\s*MOL WT")
        gamma = first_column(r"GAMMAs")
        cstar = first_column(r"CSTAR,\s*M/SEC")

        found = {k: v for k, v in
                 dict(p_bar=p_bar, T=T, MW=MW, gamma=gamma, cstar=cstar).items()
                 if v is not None}
        if None in (p_bar, T, MW, gamma):
            raise ValueError(
                f"{path}: could not read a CEA rocket case. Found {found}. "
                f"Expected lines for 'P, BAR', 'T, K', 'M, (1/n)' and "
                f"'GAMMAs'. Use a CSV table instead if this file is an "
                f"unusual CEA variant.")

        comp = {}
        m = re.search(r"MOLE FRACTIONS(.*?)(?:\n\s*\n|\Z)", text, re.S)
        if m:
            for line in m.group(1).splitlines():
                mm = re.match(r"\s*\*?([A-Za-z0-9()*+\-]+)\s+([\d.Ee+-]+)\s*$", line)
                if mm:
                    comp[mm.group(1)] = float(mm.group(2))
        return cls(name or path.stem, p=p_bar * 1e5, T=T, MW=MW, gamma=gamma,
                   cstar=cstar, composition=comp, source=str(path))

    # -- helpers ----------------------------------------------------------
    def summary(self, p0=None):
        c = self.chamber_conditions(p0)
        L = [f"Propellant: {self.name}"]
        if self.source:
            L.append(f"  source            {self.source}")
        L.append(f"  chamber pressure  {c['p'] / 1e5:10.4f}  bar"
                 + ("" if self.p.size > 1 else "   (pressure-independent table)"))
        L.append(f"  temperature       {c['T']:10.1f}  K")
        L.append(f"  molecular weight  {c['MW']:10.3f}  kg/kmol")
        L.append(f"  gamma             {c['gamma']:10.4f}")
        L.append(f"  c*                {c['cstar']:10.1f}  m/s"
                 + ("" if self.cstar_table is not None else "   (computed, not from CEA)"))
        if self.composition:
            top = sorted(self.composition.items(), key=lambda kv: -kv[1])[:6]
            L.append("  major species     "
                     + ", ".join(f"{k} {v:.3f}" for k, v in top))
        return "\n".join(L)

    def __repr__(self):
        return (f"TabulatedPropellant({self.name!r}, "
                f"{self.p.size} pressure point(s))")


# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------
def available(directory=None):
    """Names of propellant tables found on disk."""
    d = Path(directory or PROPELLANT_DIR)
    return sorted(p.stem for p in d.glob("*.csv")) if d.is_dir() else []


def load(name, directory=None):
    """Load a propellant by name from the propellants directory."""
    d = Path(directory or PROPELLANT_DIR)
    path = d / f"{name}.csv"
    if not path.exists():
        have = available(d)
        raise FileNotFoundError(
            f"no propellant table {name!r} in {d}. "
            + (f"Available: {have}." if have else
               "The directory is empty -- add a CSV, or call "
               "write_template() to start one."))
    return TabulatedPropellant.from_csv(path, name=name)


TEMPLATE = """\
# Chamber conditions for {name}.
#
# One row makes the propellant pressure-independent; add rows to interpolate.
# Fill these from a CEA rocket run at your chamber conditions -- do not guess.
# cstar_m_s is optional but preferred: CEA accounts for shifting equilibrium
# through the throat, which the closed-form expression does not.
# mu_Pa_s and Pr are optional; omit them to use a power law anchored at T_K.
# x_* columns are optional and carried for reference only.
#
p_bar,T_K,MW,gamma,cstar_m_s,mu_Pa_s,Pr
,,,,,,
"""


def write_template(name, directory=None):
    """Write a blank propellant table to fill in from CEA."""
    d = Path(directory or PROPELLANT_DIR)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.csv"
    if path.exists():
        raise FileExistsError(f"{path} already exists")
    path.write_text(TEMPLATE.format(name=name), encoding="utf-8")
    return path


def from_mechanism(mech, X, name="hydrazine", p_bar=None):
    """Build a table from a reacting mechanism, for cross-checking.

    Lets a modelled propellant and a tabulated one be compared through exactly
    the same interface -- useful for confirming the tabulated path reproduces
    what the hydrazine mechanism already computes.
    """
    c = mech.chamber_conditions(X)
    return TabulatedPropellant(
        name, p=(p_bar or 1.0) * 1e5, T=c["T"], MW=c["MW"], gamma=c["gamma"],
        cstar=c["cstar"], source=f"{type(mech).__name__}, X={X}")


#: Clearer names at package level.
load_propellant = load
available_propellants = available
