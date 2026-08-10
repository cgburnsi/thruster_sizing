"""Structured, body-fitted axisymmetric grid generation and metrics.

Index conventions
-----------------
* ``i`` runs streamwise (axial), ``j`` runs radially outward.
* ``j = 0`` lies on the symmetry axis (r = 0); ``j = nj`` lies on the wall.
* Nodes are ``(ni+1, nj+1)``; cells are ``(ni, nj)``.
* i-faces ``(ni+1, nj)`` separate cells in i; their unit normals point in +i.
* j-faces ``(ni, nj+1)`` separate cells in j; their unit normals point in +j.

All axisymmetric areas and volumes drop the common factor of 2*pi, i.e.
``S = integral(r dl)`` and ``V = integral(r dA)``. The factor cancels
everywhere in the finite-volume update, and is reinstated explicitly when
integrating physical quantities such as mass flow and thrust.
"""
import numpy as np


def _equidistribute(x_lo, x_hi, n, weight, n_probe=20001):
    """Place ``n+1`` nodes on [x_lo, x_hi] so that ``weight(x) dx`` is equal
    between consecutive nodes."""
    xs = np.linspace(x_lo, x_hi, n_probe)
    w = np.asarray(weight(xs), dtype=float)
    w = np.maximum(w, 1e-12)
    # cumulative trapezoidal integral of the weight
    c = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(xs))])
    c /= c[-1]
    return np.interp(np.linspace(0.0, 1.0, n + 1), c, xs)


def _wall_cluster(n, beta):
    """Normalised radial distribution on [0, 1], clustered toward 1 (wall)."""
    if beta <= 1e-6:
        return np.linspace(0.0, 1.0, n + 1)
    s = np.linspace(0.0, 1.0, n + 1)
    eta = np.tanh(beta * s) / np.tanh(beta)
    eta[0], eta[-1] = 0.0, 1.0
    return eta


def _solve_beta(n, first_spacing):
    """Find the tanh stretching parameter giving ``first_spacing`` (as a
    fraction of the local radius) in the cell adjacent to the wall."""
    target = float(first_spacing)
    lo, hi = 1e-4, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        eta = _wall_cluster(n, mid)
        d = eta[-1] - eta[-2]
        if d > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class Grid:
    """Holds node coordinates and all derived finite-volume metrics."""

    def __init__(self, x_nodes, r_nodes):
        self.x_n = np.asarray(x_nodes, dtype=float)
        self.r_n = np.asarray(r_nodes, dtype=float)
        if self.x_n.shape != self.r_n.shape or self.x_n.ndim != 2:
            raise ValueError("node arrays must be 2-D and the same shape")
        self.ni = self.x_n.shape[0] - 1
        self.nj = self.x_n.shape[1] - 1
        self._compute_metrics()

    # -- construction -----------------------------------------------------
    @classmethod
    def from_contour(cls, contour, ni, nj, wall_spacing=0.01,
                     w_slope=2.0, w_throat=3.0):
        """Algebraic grid: vertical i-lines, radial lines clustered at the wall.

        Parameters
        ----------
        contour : NozzleContour
        ni, nj : int
            Cell counts in the axial and radial directions.
        wall_spacing : float
            Height of the first cell off the wall as a fraction of the local
            radius. 0.01 with nj=80 gives a strongly stretched, boundary-layer
            resolving distribution.
        w_slope, w_throat : float
            Axial clustering weights. ``w_slope`` adds points where the wall
            turns; ``w_throat`` adds points where the radius is small.
        """
        Rt = contour.r_throat

        def weight(x):
            r = contour.r_wall(x)
            h = 1e-6 * (contour.x_end - contour.x_start)
            drdx = (contour.r_wall(x + h) - contour.r_wall(x - h)) / (2 * h)
            return 1.0 + w_slope * np.abs(drdx) + w_throat * (Rt / r)

        xs = _equidistribute(contour.x_start, contour.x_end, ni, weight)
        beta = _solve_beta(nj, wall_spacing)
        eta = _wall_cluster(nj, beta)

        rw = contour.r_wall(xs)                      # (ni+1,)
        x_nodes = np.repeat(xs[:, None], nj + 1, axis=1)
        r_nodes = rw[:, None] * eta[None, :]
        g = cls(x_nodes, r_nodes)
        g.contour = contour
        g.beta = beta
        return g

    # -- metrics ----------------------------------------------------------
    def _compute_metrics(self):
        x, r = self.x_n, self.r_n

        # Cell vertices, counter-clockwise in the (x, r) plane
        x0, r0 = x[:-1, :-1], r[:-1, :-1]
        x1, r1 = x[1:, :-1], r[1:, :-1]
        x2, r2 = x[1:, 1:], r[1:, 1:]
        x3, r3 = x[:-1, 1:], r[:-1, 1:]

        cx = [x0, x1, x2, x3]
        cr = [r0, r1, r2, r3]

        cross = np.zeros_like(x0)
        cxi = np.zeros_like(x0)
        cri = np.zeros_like(x0)
        vol = np.zeros_like(x0)
        for k in range(4):
            a, b = k, (k + 1) % 4
            cr_k = cx[a] * cr[b] - cx[b] * cr[a]
            cross += cr_k
            cxi += (cx[a] + cx[b]) * cr_k
            cri += (cr[a] + cr[b]) * cr_k
            vol += (cx[b] - cx[a]) * (cr[a] ** 2 + cr[a] * cr[b] + cr[b] ** 2) / 3.0

        self.A_planar = 0.5 * cross                 # planar cell area [m^2]
        self.V = -0.5 * vol                         # integral(r dA)  [m^3 / 2pi]
        self.xc = cxi / (6.0 * self.A_planar)
        self.rc = cri / (6.0 * self.A_planar)

        if np.any(self.A_planar <= 0.0) or np.any(self.V <= 0.0):
            raise ValueError("grid contains inverted or zero-volume cells")

        # i-faces: node (i,j) -> (i,j+1); normal points in +i
        dxi = x[:, 1:] - x[:, :-1]
        dri = r[:, 1:] - r[:, :-1]
        Li = np.hypot(dxi, dri)
        self.L_i = Li
        self.nx_i = dri / Li
        self.nr_i = -dxi / Li
        self.S_i = Li * 0.5 * (r[:, 1:] + r[:, :-1])
        self.xf_i = 0.5 * (x[:, 1:] + x[:, :-1])
        self.rf_i = 0.5 * (r[:, 1:] + r[:, :-1])

        # j-faces: node (i,j) -> (i+1,j); normal points in +j
        dxj = x[1:, :] - x[:-1, :]
        drj = r[1:, :] - r[:-1, :]
        Lj = np.hypot(dxj, drj)
        self.L_j = Lj
        self.nx_j = -drj / Lj
        self.nr_j = dxj / Lj
        self.S_j = Lj * 0.5 * (r[1:, :] + r[:-1, :])
        self.xf_j = 0.5 * (x[1:, :] + x[:-1, :])
        self.rf_j = 0.5 * (r[1:, :] + r[:-1, :])

        # Directional face-area vectors used by the local time-step estimate
        self.Sx_I = 0.5 * (self.nx_i[:-1] * self.S_i[:-1] + self.nx_i[1:] * self.S_i[1:])
        self.Sr_I = 0.5 * (self.nr_i[:-1] * self.S_i[:-1] + self.nr_i[1:] * self.S_i[1:])
        self.Sx_J = 0.5 * (self.nx_j[:, :-1] * self.S_j[:, :-1] + self.nx_j[:, 1:] * self.S_j[:, 1:])
        self.Sr_J = 0.5 * (self.nr_j[:, :-1] * self.S_j[:, :-1] + self.nr_j[:, 1:] * self.S_j[:, 1:])

    # -- diagnostics ------------------------------------------------------
    def wall_normal_spacing(self):
        """Physical height of the cell layer adjacent to the wall [m]."""
        return np.hypot(self.xc[:, -1] - self.xf_j[:, -1],
                        self.rc[:, -1] - self.rf_j[:, -1]) * 2.0

    def summary(self):
        h = self.wall_normal_spacing()
        dx = np.diff(self.x_n[:, 0])
        ar = self.A_planar / np.minimum(self.L_i[:-1], self.L_i[1:]) ** 2
        return (
            f"Grid {self.ni} x {self.nj} = {self.ni * self.nj} cells\n"
            f"  axial spacing      : {dx.min():.3e} .. {dx.max():.3e} m\n"
            f"  wall cell height   : {h.min():.3e} .. {h.max():.3e} m\n"
            f"  cell volume (/2pi) : {self.V.min():.3e} .. {self.V.max():.3e} m^3\n"
            f"  max aspect ratio   : {ar.max():.1f}"
        )

    def cell_centre_delta_i(self):
        """Vector between the cell centres straddling each i-face (interior
        faces only), returned as (dx, dr, dist)."""
        dx = self.xc[1:, :] - self.xc[:-1, :]
        dr = self.rc[1:, :] - self.rc[:-1, :]
        return dx, dr, np.hypot(dx, dr)

    def cell_centre_delta_j(self):
        dx = self.xc[:, 1:] - self.xc[:, :-1]
        dr = self.rc[:, 1:] - self.rc[:, :-1]
        return dx, dr, np.hypot(dx, dr)
