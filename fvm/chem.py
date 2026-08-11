"""Multi-species thermodynamics for reacting flow.

The nozzle solver treats the gas as a single calorically perfect fluid, which
is fine downstream of a catalyst bed where composition is frozen. Inside the
bed it is not: hydrazine decomposition changes the mixture continuously, and
the ammonia dissociation fraction sets T0, molecular weight and gamma
together. That needs real species thermodynamics.

Properties come from NASA 7-coefficient polynomials, the same form used by
CHEMKIN and by NASA CEA's thermo data::

    cp/Ru = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
    h/(Ru T) = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
    s/Ru = a1 lnT + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7

with separate coefficient sets below and above a common temperature.

On the built-in data
--------------------
The table here covers the species hydrazine decomposition needs. Coefficients
are transcribed from standard thermodynamic databases, and **transcription is
exactly the kind of thing that fails silently** -- a wrong digit in a5 shifts
cp by a few percent at high temperature and nothing crashes. Every species is
therefore pinned by tests against independently known cp values at several
temperatures, so a bad coefficient fails the suite rather than quietly
corrupting a bed calculation.

For authoritative work, load your own data with :func:`load_chemkin` rather
than trusting the built-in table.
"""
import numpy as np

RU = 8314.462618          # J/(kmol*K)


class Species:
    """One chemical species with NASA-7 polynomial thermodynamics.

    Parameters
    ----------
    name : str
    MW : float
        Molecular weight [kg/kmol].
    low, high : sequence of 7 floats
        NASA-7 coefficients below and above ``T_mid``.
    T_range : (float, float, float)
        (T_min, T_mid, T_max) in K.
    """

    def __init__(self, name, MW, low, high, T_range=(200.0, 1000.0, 6000.0)):
        self.name = name
        self.MW = float(MW)
        self.low = np.asarray(low, dtype=float)
        self.high = np.asarray(high, dtype=float)
        if self.low.size != 7 or self.high.size != 7:
            raise ValueError("NASA-7 needs exactly 7 coefficients per range")
        self.T_min, self.T_mid, self.T_max = map(float, T_range)
        self.R = RU / self.MW          # specific gas constant [J/(kg*K)]

    def _coeffs(self, T):
        T = np.asarray(T, dtype=float)
        lo = self.low.reshape(7, *([1] * T.ndim))
        hi = self.high.reshape(7, *([1] * T.ndim))
        return np.where(T < self.T_mid, lo, hi)

    def cp_mole(self, T):
        """Molar heat capacity at constant pressure [J/(kmol*K)]."""
        T = np.asarray(T, dtype=float)
        a = self._coeffs(T)
        return RU * (a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] * T ** 4)

    def h_mole(self, T):
        """Molar enthalpy [J/kmol], including enthalpy of formation."""
        T = np.asarray(T, dtype=float)
        a = self._coeffs(T)
        return RU * T * (a[0] + a[1] * T / 2.0 + a[2] * T ** 2 / 3.0
                         + a[3] * T ** 3 / 4.0 + a[4] * T ** 4 / 5.0 + a[5] / T)

    def s_mole(self, T):
        """Standard-state molar entropy [J/(kmol*K)] at 1 bar."""
        T = np.asarray(T, dtype=float)
        a = self._coeffs(T)
        return RU * (a[0] * np.log(T) + a[1] * T + a[2] * T ** 2 / 2.0
                     + a[3] * T ** 3 / 3.0 + a[4] * T ** 4 / 4.0 + a[6])

    # mass-specific forms
    def cp(self, T):
        """Specific heat at constant pressure [J/(kg*K)]."""
        return self.cp_mole(T) / self.MW

    def h(self, T):
        """Specific enthalpy [J/kg], including enthalpy of formation."""
        return self.h_mole(T) / self.MW

    def h_formation(self):
        """Standard enthalpy of formation at 298.15 K [J/kmol]."""
        return self.h_mole(298.15)

    def __repr__(self):
        return f"Species({self.name!r}, MW={self.MW:.4f})"


class Mixture:
    """An ideal-gas mixture held as mass fractions over a fixed species list."""

    def __init__(self, species, Y=None):
        self.species = list(species)
        self.names = [s.name for s in self.species]
        self.MW_k = np.array([s.MW for s in self.species])
        self.n = len(self.species)
        self.Y = np.zeros(self.n) if Y is None else self.normalized(Y)

    # -- composition helpers ----------------------------------------------
    def index(self, name):
        return self.names.index(name)

    @staticmethod
    def normalized(Y):
        Y = np.asarray(Y, dtype=float)
        tot = Y.sum(axis=0)
        return Y / np.where(np.abs(tot) < 1e-300, 1.0, tot)

    def from_moles(self, moles):
        """Mass fractions from a dict or array of mole amounts."""
        if isinstance(moles, dict):
            x = np.zeros(self.n)
            for k, v in moles.items():
                x[self.index(k)] = v
        else:
            x = np.asarray(moles, dtype=float)
        return self.normalized(x * self.MW_k)

    def mole_fractions(self, Y=None):
        Y = self.Y if Y is None else Y
        x = np.asarray(Y, dtype=float) / self.MW_k.reshape(-1, *([1] * (np.ndim(Y) - 1)))
        return self.normalized(x)

    # -- mixture properties -----------------------------------------------
    def MW(self, Y=None):
        """Mixture molecular weight [kg/kmol]."""
        Y = self.Y if Y is None else Y
        shape = (-1,) + (1,) * (np.ndim(Y) - 1)
        return 1.0 / np.sum(np.asarray(Y) / self.MW_k.reshape(shape), axis=0)

    def R(self, Y=None):
        """Specific gas constant [J/(kg*K)]."""
        return RU / self.MW(Y)

    def cp(self, T, Y=None):
        """Mixture specific heat [J/(kg*K)]."""
        Y = self.Y if Y is None else Y
        return sum(Y[k] * s.cp(T) for k, s in enumerate(self.species))

    def h(self, T, Y=None):
        """Mixture specific enthalpy [J/kg], including formation enthalpies."""
        Y = self.Y if Y is None else Y
        return sum(Y[k] * s.h(T) for k, s in enumerate(self.species))

    def gamma(self, T, Y=None):
        cp = self.cp(T, Y)
        return cp / (cp - self.R(Y))

    def sound_speed(self, T, Y=None):
        return np.sqrt(self.gamma(T, Y) * self.R(Y) * T)

    def temperature_from_h(self, h_target, Y=None, T_guess=1000.0,
                           tol=1e-8, max_iter=100):
        """Invert h(T) for temperature by Newton iteration.

        Needed because the energy equation carries enthalpy while the rate
        expressions want temperature; with variable cp there is no closed form.
        """
        Y = self.Y if Y is None else Y
        T = np.full_like(np.asarray(h_target, dtype=float), float(T_guess))
        for _ in range(max_iter):
            f = self.h(T, Y) - h_target
            dfdT = np.maximum(self.cp(T, Y), 1e-6)
            dT = -f / dfdT
            dT = np.clip(dT, -500.0, 500.0)          # keep early steps sane
            T = np.clip(T + dT, 50.0, 6000.0)
            if np.all(np.abs(dT) < tol * np.maximum(T, 1.0)):
                return T
        return T

    def __repr__(self):
        return f"Mixture({', '.join(self.names)})"


_ATOMIC = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994,
           "AR": 39.948, "HE": 4.002602, "S": 32.065}


# ---------------------------------------------------------------------------
# Built-in species data. Pinned by tests -- see the module docstring.
# ---------------------------------------------------------------------------
#: Species are stored with their *elemental composition*, not a molecular
#: weight, and MW is computed from a single atomic table. Quoting MWs directly
#: mixes atomic-weight conventions between species -- an early revision had
#: NH3 on H = 1.00797 and N2H4 on H = 1.00794, which broke mass conservation
#: in the reaction stoichiometry by a few parts per million. Deriving them
#: makes sum(nu_i MW_i) = 0 exact by construction.
_DATA = {
    "N2": ({"N": 2},
           [3.298677e0, 1.4082404e-3, -3.963222e-6, 5.641515e-9, -2.444854e-12,
            -1.0208999e3, 3.950372e0],
           [2.92664e0, 1.4879768e-3, -5.68476e-7, 1.0097038e-10, -6.753351e-15,
            -9.227977e2, 5.980528e0],
           (300.0, 1000.0, 5000.0)),
    "H2": ({"H": 2},
           [3.298124e0, 8.249442e-4, -8.143015e-7, -9.475434e-11, 4.134872e-13,
            -1.012521e3, -3.294094e0],
           [2.991423e0, 7.000644e-4, -5.633829e-8, -9.231578e-12, 1.582752e-15,
            -8.35034e2, -1.35511e0],
           (300.0, 1000.0, 5000.0)),
    "NH3": ({"N": 1, "H": 3},
            [4.2860274e0, -4.660523e-3, 2.1715133e-5, -2.2808887e-8, 8.2638046e-12,
             -6.7417285e3, -6.2537277e-1],
            [2.7170969e0, 5.5685644e-3, -1.7688659e-6, 2.6741782e-10, -1.5273113e-14,
             -6.5845128e3, 6.0928908e0],
            (200.0, 1000.0, 6000.0)),
    "H2O": ({"H": 2, "O": 1},
            [4.19864056e0, -2.0364341e-3, 6.52040211e-6, -5.48797062e-9, 1.77197817e-12,
             -3.02937267e4, -8.49032208e-1],
            [3.03399249e0, 2.17691804e-3, -1.64072518e-7, -9.7041987e-11, 1.68200992e-14,
             -3.00042971e4, 4.9667701e0],
            (200.0, 1000.0, 3500.0)),
    "O2": ({"O": 2},
           [3.78245636e0, -2.99673416e-3, 9.84730201e-6, -9.68129509e-9, 3.24372837e-12,
            -1.06394356e3, 3.65767573e0],
           [3.28253784e0, 1.48308754e-3, -7.57966669e-7, 2.09470555e-10, -2.16717794e-14,
            -1.08845772e3, 5.45323129e0],
           (200.0, 1000.0, 3500.0)),
}

#: Gas-phase hydrazine. Coefficients are less standardised than the product
#: species; the enthalpy of formation is the number that matters most for bed
#: energetics, and it is pinned by test against the accepted +95.4 kJ/mol.
_DATA["N2H4"] = (
    {"N": 2, "H": 4},
    [3.83472149e0, -6.49129555e-4, 3.76848463e-5, -5.00709182e-8, 2.03362064e-11,
     1.00893925e4, 5.7527203e0],
    [4.93957357e0, 8.75017187e-3, -2.99399058e-6, 4.67278418e-10, -2.73068599e-14,
     9.28264683e3, 2.6944021e0],
    (200.0, 1000.0, 6000.0))


def composition(name):
    """Elemental composition of a built-in species, e.g. ``{"N": 2, "H": 4}``."""
    try:
        return dict(_DATA[name][0])
    except KeyError:
        raise KeyError(f"no built-in data for {name!r}; "
                       f"have {sorted(_DATA)} (or use load_chemkin)") from None


def molecular_weight(comp):
    """Molecular weight [kg/kmol] from an elemental composition."""
    return sum(_ATOMIC[el.upper()] * n for el, n in comp.items())


def species(name):
    """Look up a built-in species by name."""
    try:
        comp, low, high, rng = _DATA[name]
    except KeyError:
        raise KeyError(f"no built-in data for {name!r}; "
                       f"have {sorted(_DATA)} (or use load_chemkin)") from None
    sp = Species(name, molecular_weight(comp), low, high, rng)
    sp.composition = dict(comp)
    return sp


def mixture(names, Y=None):
    """Build a :class:`Mixture` from built-in species names."""
    return Mixture([species(n) for n in names], Y)


def load_chemkin(path, names=None):
    """Read species from a CHEMKIN-format thermo file (NASA-7, 4 lines each).

    Use this in preference to the built-in table when the numbers matter --
    CEA's ``thermo.inp`` and Burcat's tables are both available in this format.
    """
    out = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()
                 and not ln.lstrip().startswith("!")]
    i = 0
    while i + 3 < len(lines):
        l1 = lines[i]
        if len(l1) < 80 or l1[79] != "1":
            i += 1
            continue
        name = l1[:18].split()[0]
        try:
            tlo, thi, tmid = (float(x) for x in l1[45:73].split()[:3])
            nums = []
            for ln in lines[i + 1:i + 4]:
                nums += [float(ln[j:j + 15]) for j in range(0, 75, 15)
                         if ln[j:j + 15].strip()]
            MW = _elemental_MW(l1[24:44])
            out[name] = Species(name, MW, nums[7:14], nums[0:7], (tlo, tmid, thi))
        except (ValueError, IndexError):
            pass
        i += 4
    if names is not None:
        missing = [n for n in names if n not in out]
        if missing:
            raise KeyError(f"{path} lacks species {missing}")
        return [out[n] for n in names]
    return out


def _elemental_MW(field):
    """Molecular weight from the CHEMKIN element field (4 x 'EE nnn')."""
    MW = 0.0
    for j in range(0, 20, 5):
        chunk = field[j:j + 5]
        el = chunk[:2].strip().upper()
        cnt = chunk[2:].strip()
        if el and cnt:
            MW += _ATOMIC.get(el, 0.0) * float(cnt)
    return MW
