"""Parametric nozzle contours.

A contour is defined by its wall radius as a function of axial station,
``r_wall(x)``, on ``[x_start, x_end]`` with the throat at ``x = 0``.

Conventions
-----------
* ``x`` increases in the flow direction; the throat sits at ``x = 0``.
* The chamber occupies ``x < 0``, the divergent section ``x > 0``.
* All lengths are in metres.
"""
import numpy as np


class NozzleContour:
    """Base class: holds the assembled piecewise contour as breakpoints.

    Subclasses populate ``self._x`` and ``self._r`` with a densely sampled
    polyline; ``r_wall`` then interpolates it. Sampling densely (rather than
    evaluating the analytic pieces) keeps the interface uniform and lets
    measured/CAD contours be dropped in through :meth:`from_points`.
    """

    def __init__(self, x, r):
        x = np.asarray(x, dtype=float)
        r = np.asarray(r, dtype=float)
        order = np.argsort(x)
        self._x = x[order]
        self._r = r[order]
        if np.any(self._r <= 0.0):
            raise ValueError("contour radius must be strictly positive")

    # -- geometry queries -------------------------------------------------
    @property
    def x_start(self):
        return float(self._x[0])

    @property
    def x_end(self):
        return float(self._x[-1])

    @property
    def r_throat(self):
        return float(self._r.min())

    @property
    def x_throat(self):
        return float(self._x[np.argmin(self._r)])

    @property
    def r_exit(self):
        return float(self._r[-1])

    @property
    def r_chamber(self):
        return float(self._r[0])

    @property
    def area_ratio(self):
        return (self.r_exit / self.r_throat) ** 2

    @property
    def contraction_ratio(self):
        return (self.r_chamber / self.r_throat) ** 2

    def r_wall(self, x):
        return np.interp(x, self._x, self._r)

    def area(self, x):
        return np.pi * self.r_wall(x) ** 2

    def polyline(self):
        return self._x.copy(), self._r.copy()

    @classmethod
    def from_points(cls, x, r):
        """Build a contour from a measured or CAD-exported polyline."""
        return cls(x, r)

    def __repr__(self):
        return (f"{type(self).__name__}(L={self.x_end - self.x_start:.4e} m, "
                f"Rt={self.r_throat:.4e} m, eps={self.area_ratio:.1f}, "
                f"CR={self.contraction_ratio:.1f})")


class ConicalNozzle(NozzleContour):
    """Chamber + convergent cone + throat arcs + divergent cone.

    Parameters
    ----------
    r_chamber, r_throat : float
        Chamber and throat radii [m].
    area_ratio : float
        Exit-to-throat area ratio, A_e/A_t.
    L_chamber : float
        Length of the constant-radius chamber section [m].
    theta_conv, theta_div : float
        Convergent and divergent cone half-angles [deg].
    Ru_frac, Rd_frac : float
        Throat blend-arc radii as multiples of the throat radius, upstream
        and downstream respectively. Typical values are 1.5 and 0.4.
    n_sample : int
        Number of points used to discretise the contour.
    """

    def __init__(self, r_chamber, r_throat, area_ratio, L_chamber,
                 theta_conv=35.0, theta_div=15.0,
                 Ru_frac=1.5, Rd_frac=0.4, n_sample=4001):
        Rt = float(r_throat)
        Rc = float(r_chamber)
        Re = Rt * np.sqrt(area_ratio)
        if Rc <= Rt:
            raise ValueError("chamber radius must exceed throat radius")
        if Re <= Rt:
            raise ValueError("area_ratio must exceed 1")

        tc = np.radians(theta_conv)
        td = np.radians(theta_div)
        Ru = Ru_frac * Rt
        Rd = Rd_frac * Rt

        # Upstream throat arc, centre at (0, Rt + Ru); phi measured from throat
        x1 = -Ru * np.sin(tc)
        r1 = Rt + Ru * (1.0 - np.cos(tc))
        # Convergent cone runs from the chamber radius down to (x1, r1)
        x0 = x1 - (Rc - r1) / np.tan(tc)
        if x0 >= x1:
            raise ValueError("convergent cone degenerate: increase r_chamber "
                             "or reduce Ru_frac/theta_conv")
        x_in = x0 - L_chamber

        # Downstream throat arc, centre at (0, Rt + Rd)
        x2 = Rd * np.sin(td)
        r2 = Rt + Rd * (1.0 - np.cos(td))
        # Divergent cone out to the exit radius
        x3 = x2 + (Re - r2) / np.tan(td)

        self.theta_conv = theta_conv
        self.theta_div = theta_div
        self.x_inlet, self.x_conv_start = x_in, x0
        self.x_arc_up, self.x_arc_dn = x1, x2
        self.x_exit_plane = x3

        # Force x = 0 and both arc tangency points into the sample set, so
        # that r_throat is the true throat radius and not whatever the nearest
        # sample happened to land on.
        x = np.union1d(np.linspace(x_in, x3, int(n_sample)),
                       np.array([x0, x1, 0.0, x2]))
        r = np.empty_like(x)

        m0 = x <= x0                                # chamber
        m1 = (x > x0) & (x <= x1)                   # convergent cone
        m2 = (x > x1) & (x <= 0.0)                  # upstream throat arc
        m3 = (x > 0.0) & (x <= x2)                  # downstream throat arc
        m4 = x > x2                                 # divergent cone

        r[m0] = Rc
        r[m1] = r1 + (x1 - x[m1]) * np.tan(tc)
        r[m2] = Rt + Ru - np.sqrt(np.maximum(Ru * Ru - x[m2] ** 2, 0.0))
        r[m3] = Rt + Rd - np.sqrt(np.maximum(Rd * Rd - x[m3] ** 2, 0.0))
        r[m4] = r2 + (x[m4] - x2) * np.tan(td)

        super().__init__(x, r)


class BellNozzle(NozzleContour):
    """Rao-style bell: same convergent side, quadratic-Bezier divergent side.

    The divergent wall leaves the downstream throat arc at ``theta_n`` and
    arrives at the exit plane at ``theta_e``. Its length is ``bell_frac``
    times the length of a 15-degree cone of the same area ratio.
    """

    def __init__(self, r_chamber, r_throat, area_ratio, L_chamber,
                 theta_conv=35.0, theta_n=30.0, theta_e=8.0,
                 bell_frac=0.8, Ru_frac=1.5, Rd_frac=0.382, n_sample=4001):
        Rt = float(r_throat)
        Rc = float(r_chamber)
        Re = Rt * np.sqrt(area_ratio)

        tc = np.radians(theta_conv)
        tn = np.radians(theta_n)
        te = np.radians(theta_e)
        Ru = Ru_frac * Rt
        Rd = Rd_frac * Rt

        x1 = -Ru * np.sin(tc)
        r1 = Rt + Ru * (1.0 - np.cos(tc))
        x0 = x1 - (Rc - r1) / np.tan(tc)
        x_in = x0 - L_chamber

        # Bezier start point: end of the downstream throat arc
        xN = Rd * np.sin(tn)
        rN = Rt + Rd * (1.0 - np.cos(tn))
        # Exit point, placed at bell_frac of the equivalent 15-deg cone length
        L_cone15 = (Re - Rt) / np.tan(np.radians(15.0))
        xE = bell_frac * L_cone15
        if xE <= xN:
            raise ValueError("bell_frac too small: divergent section degenerate")
        rE = Re

        # Control point = intersection of the two tangent lines
        # line N: r = rN + tan(tn)*(x - xN);  line E: r = rE + tan(te)*(x - xE)
        denom = np.tan(tn) - np.tan(te)
        if abs(denom) < 1e-12:
            raise ValueError("theta_n and theta_e must differ")
        xC = (rE - rN + np.tan(tn) * xN - np.tan(te) * xE) / denom
        rC = rN + np.tan(tn) * (xC - xN)

        self.theta_conv = theta_conv
        self.theta_n, self.theta_e = theta_n, theta_e
        self.x_inlet, self.x_conv_start = x_in, x0
        self.x_arc_up, self.x_arc_dn = x1, xN
        self.x_exit_plane = xE

        # Sample: convergent side on a uniform x grid, bell on the Bezier
        n_up = int(n_sample) // 2
        xu = np.union1d(np.linspace(x_in, xN, n_up), np.array([x0, x1, 0.0]))
        ru = np.empty_like(xu)
        m0 = xu <= x0
        m1 = (xu > x0) & (xu <= x1)
        m2 = (xu > x1) & (xu <= 0.0)
        m3 = xu > 0.0
        ru[m0] = Rc
        ru[m1] = r1 + (x1 - xu[m1]) * np.tan(tc)
        ru[m2] = Rt + Ru - np.sqrt(np.maximum(Ru * Ru - xu[m2] ** 2, 0.0))
        ru[m3] = Rt + Rd - np.sqrt(np.maximum(Rd * Rd - xu[m3] ** 2, 0.0))

        t = np.linspace(0.0, 1.0, int(n_sample) - n_up)[1:]
        xb = (1 - t) ** 2 * xN + 2 * (1 - t) * t * xC + t ** 2 * xE
        rb = (1 - t) ** 2 * rN + 2 * (1 - t) * t * rC + t ** 2 * rE

        super().__init__(np.concatenate([xu, xb]), np.concatenate([ru, rb]))
