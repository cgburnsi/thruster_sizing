"""Fit catalyst bed kinetics to measured data.

Why this exists
---------------
The bed model's structure is well founded; its rate constants are not. Kesten
fitted them to engine data in 1967 and was explicit about their standing -- the
hydrazine activation energy was "chosen rather arbitrarily", the hydrogen
inhibition order was measured on platinum and assumed to transfer to Shell 405,
and the ammonia pre-exponential carries a factor-of-three band. Between
``A_NH3`` and ``n_H2`` the predicted dissociation spans more than an order of
magnitude.

So the honest use of this model is to fit those two numbers to *your* bed and
then predict with them. This module does the fitting; it holds no GUI code, so
it can be tested and scripted on its own.

What can be fitted
------------------
``A_NH3``    ammonia pre-exponential, in Kesten's units. The single parameter
             with the most leverage over predicted dissociation.
``n_H2``     hydrogen inhibition order. Unresolved in Kesten's own source --
             1.0 in ``MAIN.f``, 1.6 in ``PARAM.f`` and his published Eq. 43.
``T_vapor``  where the vapour region begins. A modelling input rather than a
             measurement, and the rates are sensitive to it.
``A_N2H4``   hydrazine pre-exponential. Included for completeness, but fitting
             it is usually pointless and that is worth seeing: the reaction is
             diffusion-controlled, so the objective is nearly flat in it.

Sensitivity is reported alongside any fit, because a parameter the data cannot
constrain will otherwise come back with a confident-looking value.
"""
import csv
from pathlib import Path

import numpy as np

from .catbed import CatalystBed
from .mechanism import HydrazineShell405, RANKINE_TO_KELVIN
from .plugflow import PlugFlowReactor
from .thruster import vapor_region_inlet

FT = 0.3048
PSI = 6894.757293
LB = 0.45359237

#: Parameters that may be fitted, with sane bounds and whether they are
#: explored logarithmically.
PARAMETERS = {
    "A_NH3": dict(default=1.0e11, bounds=(1e9, 1e13), log=True,
                  label="A_NH3 (Kesten units)"),
    "n_H2": dict(default=1.0, bounds=(0.0, 2.0), log=False,
                 label="H2 inhibition order"),
    "T_vapor": dict(default=455.6, bounds=(380.0, 900.0), log=False,
                    label="Vapour region inlet [K]"),
    "A_N2H4": dict(default=1.0e10, bounds=(1e8, 1e13), log=True,
                   label="A_N2H4 (Kesten units)"),
}


class BedData:
    """Measured profiles along a catalyst bed.

    Expected columns, matched case-insensitively with flexible names:

    ==================  ===========================================
    position            ``z``, ``z_mm``, ``z_m``, ``z_ft``, ``x``
    gas temperature     ``T``, ``T_gas``, ``T_K``, ``T_degR``, ``T_degC``
    dissociation        ``X``, ``dissociation``, ``frac3d``
    pressure            ``p``, ``p_bar``, ``p_psia``, ``p_Pa``
    ==================  ===========================================

    Units are taken from the column suffix where present and default to SI.
    ``frac3d`` is treated as Kesten's convention and converted.
    """

    def __init__(self, z, T=None, X=None, p=None, name="measured"):
        self.z = np.asarray(z, dtype=float)
        self.T = None if T is None else np.asarray(T, dtype=float)
        self.X = None if X is None else np.asarray(X, dtype=float)
        self.p = None if p is None else np.asarray(p, dtype=float)
        self.name = name
        order = np.argsort(self.z)
        self.z = self.z[order]
        for a in ("T", "X", "p"):
            v = getattr(self, a)
            if v is not None:
                setattr(self, a, v[order])

    @property
    def channels(self):
        return [a for a in ("T", "X", "p") if getattr(self, a) is not None]

    def __len__(self):
        return self.z.size

    def __repr__(self):
        return (f"BedData({self.name!r}, {len(self)} points, "
                f"channels={self.channels})")

    # -- loading ----------------------------------------------------------
    @staticmethod
    def _column(header, options):
        low = {h.strip().lower(): h for h in header}
        for want in options:
            if want in low:
                return low[want]
        return None

    @classmethod
    def from_csv(cls, path):
        path = Path(path)
        with path.open(newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            raise ValueError(f"{path} has no data rows")
        header = list(rows[0].keys())

        def col(options):
            return cls._column(header, options)

        # "x" is deliberately not a position alias: it is the conventional
        # symbol for dissociation fraction here, and accepting it for position
        # silently misreads a perfectly reasonable data file.
        z_col = col(["z_mm", "z_m", "z_ft", "z", "x_mm", "x_m", "position"])
        if z_col is None:
            raise ValueError(
                f"{path}: no position column found. Expected one of "
                f"z_mm, z_m, z_ft, z, x. Columns present: {header}")
        scale = {"z_mm": 1e-3, "z_ft": FT, "x_mm": 1e-3}.get(z_col.lower(), 1.0)  # noqa: E501
        z = np.array([float(r[z_col]) for r in rows]) * scale

        T = X = p = None
        t_col = col(["t_k", "t_gas_k", "t_gas", "temp_k", "t", "temperature",
                     "t_degr", "temp_degr", "t_degc"])
        if t_col is not None:
            raw = np.array([float(r[t_col]) for r in rows])
            key = t_col.lower()
            T = (raw * RANKINE_TO_KELVIN if "degr" in key
                 else raw + 273.15 if "degc" in key else raw)

        x_col = col(["x", "dissociation", "frac3d", "x_dissociation"])
        if x_col is not None and x_col != z_col:
            raw = np.array([float(r[x_col]) for r in rows])
            # Kesten's fraction uses a different convention; convert it.
            X = (3.0 * raw + 1.0) / 4.0 if x_col.lower() == "frac3d" else raw

        p_col = col(["p_pa", "p_bar", "p_psia", "p", "pres_psia", "pressure"])
        if p_col is not None:
            raw = np.array([float(r[p_col]) for r in rows])
            key = p_col.lower()
            p = (raw * 1e5 if "bar" in key
                 else raw * PSI if "psi" in key else raw)

        if T is None and X is None:
            raise ValueError(
                f"{path}: found a position column but nothing to fit against. "
                f"Expected a temperature or dissociation column. "
                f"Columns present: {header}")
        return cls(z, T=T, X=X, p=p, name=path.stem)


class FitCase:
    """A measured run: the bed, its operating point, and the data."""

    def __init__(self, bed, G, p_inlet, data, T_feed=298.15, Y_inlet=None,
                 z_window=None):
        self.bed = bed
        self.G = float(G)
        self.p_inlet = float(p_inlet)
        self.data = data
        self.T_feed = float(T_feed)
        self.Y_inlet = Y_inlet
        #: Axial window (z_min, z_max) of data to fit against, in metres.
        #: Excluding a region matters more than it sounds. This model begins
        #: at the vapour-region inlet and does not represent the liquid or
        #: liquid-vapour zones, so measurements there cannot be matched by any
        #: parameter value -- including them simply drags the fit off. On the
        #: Kesten demo, three of six stations lie in that zone and dominate
        #: the objective while the model reproduces the developed bed to
        #: within 1.5 K.
        self.z_window = z_window

    def mask(self):
        """Boolean mask selecting the data points inside the fitting window."""
        z = self.data.z
        if self.z_window is None:
            return np.ones(z.shape, dtype=bool)
        lo, hi = self.z_window
        lo = -np.inf if lo is None else lo
        hi = np.inf if hi is None else hi
        return (z >= lo) & (z <= hi)

    # -- the model --------------------------------------------------------
    def solve(self, params, n_output=200):
        """Integrate the bed with the given parameter values."""
        p = {k: v["default"] for k, v in PARAMETERS.items()}
        p.update(params)
        mech = HydrazineShell405(nH2=p["n_H2"], A_NH3=p["A_NH3"],
                                 A_N2H4_cat=p["A_N2H4"])
        Y0 = self.Y_inlet
        if Y0 is None:
            Y0, _ = vapor_region_inlet(mech, p["T_vapor"], self.T_feed)
        return PlugFlowReactor(mech, self.bed).solve(
            G=self.G, p_inlet=self.p_inlet, T_inlet=p["T_vapor"],
            Y_inlet=Y0, n_output=n_output)

    # -- comparison -------------------------------------------------------
    def predict_at_data(self, solution):
        """Model values interpolated onto the measurement positions."""
        out = {}
        if self.data.T is not None:
            out["T"] = np.interp(self.data.z, solution.z, solution.T_gas)
        if self.data.X is not None:
            out["X"] = np.interp(self.data.z, solution.z, solution.dissociation)
        if self.data.p is not None:
            out["p"] = np.interp(self.data.z, solution.z, solution.p)
        return out

    def _selected(self, arr):
        return np.asarray(arr)[self.mask()]

    def residuals(self, params, weights=None, solution=None):
        """Weighted residuals, one array per available channel, concatenated.

        Channels are normalised by their own spread so temperature in kelvin
        and a dimensionless dissociation fraction contribute comparably rather
        than temperature dominating by three orders of magnitude.
        """
        sol = solution if solution is not None else self.solve(params)
        pred = self.predict_at_data(sol)
        w = {"T": 1.0, "X": 1.0, "p": 0.0}
        w.update(weights or {})

        m = self.mask()
        parts = []
        for key, model in pred.items():
            meas = getattr(self.data, key)[m]
            model = model[m]
            scale = np.ptp(meas)
            if not np.isfinite(scale) or scale <= 0:
                scale = max(abs(np.mean(meas)), 1e-12)
            if w[key] > 0:
                parts.append(w[key] * (model - meas) / scale)
        if not parts:
            raise ValueError("no channels enabled for fitting")
        return np.concatenate(parts)

    def errors(self, params, solution=None):
        """Per-channel RMS error in physical units, for display."""
        sol = solution if solution is not None else self.solve(params)
        pred = self.predict_at_data(sol)
        m = self.mask()
        out = {}
        for key, model in pred.items():
            meas = getattr(self.data, key)[m]
            out[key] = float(np.sqrt(np.mean((model[m] - meas) ** 2)))
        return out, sol

    # -- built-in demo ----------------------------------------------------
    @classmethod
    def kesten_demo(cls, reference=None):
        """Kesten's own reactor output as a worked example.

        Useful because the answer is known: sweeping the parameters against
        this case recovers n_H2 = 1.0 with A_NH3 = 1e11, which is exactly what
        his MAIN.f and input deck contain. A fit that does not land there means
        the tool is wrong, not the data.
        """
        path = Path(reference) if reference else (
            Path(__file__).resolve().parents[1]
            / "docs/kesten_claude/vapor_reference.csv")
        if not path.exists():
            raise FileNotFoundError(
                f"Kesten reference data not found at {path}. It ships in "
                f"docs/kesten_claude; see docs/README.md.")
        data = BedData.from_csv(path)
        data.name = "Kesten G910461-30"
        first = data
        bed = CatalystBed.kesten_standard(diameter=0.02, length=0.25 * FT)
        # His vapour region starts already partly decomposed; use his own
        # first station rather than deriving one.
        with path.open(newline="", encoding="utf-8-sig") as fh:
            rows = sorted(csv.DictReader(fh), key=lambda r: float(r["z_ft"]))
        mech = HydrazineShell405()
        Y = np.zeros(mech.mixture.n)
        for name, key in zip(["H2", "N2", "NH3", "N2H4"], ["c1", "c2", "c3", "c4"]):
            Y[mech.mixture.index(name)] = float(rows[0][key])
        return cls(bed, G=3.0 * LB / FT ** 2,
                   p_inlet=float(rows[0]["pres_psia"]) * PSI,
                   data=first, Y_inlet=Y / Y.sum())


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------
def fit(case, fit_names=("A_NH3", "n_H2"), start=None, weights=None,
        max_nfev=60, callback=None):
    """Least-squares fit of selected parameters.

    Logarithmic parameters are optimised in log space, which matters: ``A_NH3``
    spans four decades and a linear optimiser wastes its whole budget in the
    first one.

    Returns a dict with the fitted values, the per-channel RMS errors, and a
    sensitivity estimate for each parameter -- the fractional change in
    objective for a 10% parameter change. A parameter with near-zero
    sensitivity was not constrained by the data, whatever value it came back
    with.
    """
    from scipy.optimize import least_squares

    names = list(fit_names)
    for n in names:
        if n not in PARAMETERS:
            raise ValueError(f"unknown parameter {n!r}; "
                             f"have {sorted(PARAMETERS)}")
    base = {k: v["default"] for k, v in PARAMETERS.items()}
    base.update(start or {})

    def to_x(vals):
        return np.array([np.log10(vals[n]) if PARAMETERS[n]["log"] else vals[n]
                         for n in names])

    def from_x(x):
        out = dict(base)
        for n, xi in zip(names, x):
            out[n] = 10.0 ** xi if PARAMETERS[n]["log"] else xi
        return out

    lo, hi = [], []
    for n in names:
        a, b = PARAMETERS[n]["bounds"]
        if PARAMETERS[n]["log"]:
            a, b = np.log10(a), np.log10(b)
        lo.append(a)
        hi.append(b)

    n_calls = [0]

    def residual(x):
        n_calls[0] += 1
        params = from_x(x)
        r = case.residuals(params, weights=weights)
        if callback:
            callback(n_calls[0], params, float(np.sqrt(np.mean(r ** 2))))
        return r

    res = least_squares(residual, to_x(base), bounds=(lo, hi),
                        max_nfev=max_nfev, xtol=1e-8, ftol=1e-8)
    best = from_x(res.x)
    errs, sol = case.errors(best)

    sens = {}
    r0 = float(np.sqrt(np.mean(case.residuals(best, weights=weights) ** 2)))
    for n in names:
        bumped = dict(best)
        bumped[n] = best[n] * 1.1
        a, b = PARAMETERS[n]["bounds"]
        bumped[n] = float(np.clip(bumped[n], a, b))
        r1 = float(np.sqrt(np.mean(case.residuals(bumped, weights=weights) ** 2)))
        sens[n] = abs(r1 - r0) / max(r0, 1e-30)

    return dict(params=best, errors=errs, solution=sol, rms=r0,
                sensitivity=sens, nfev=n_calls[0], success=bool(res.success))


#: Exported under a clearer name at package level.
fit_kinetics = fit
