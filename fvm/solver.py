"""Cell-centred finite-volume solver for the axisymmetric Navier-Stokes
equations.

Formulation
-----------
Working in cylindrical (x, r) with azimuthal symmetry and multiplying through
by r puts the equations in strong conservation form on the true axisymmetric
control volume::

    d(rU)/dt + d(r Fx)/dx + d(r Fr)/dr = [0, 0, p - tau_tt, 0]

so the discrete update over a cell is

    V dU/dt = -sum_faces (F_inv - F_visc) . n S + [0, 0, p - tau_tt, 0] A_planar

with ``V = integral(r dA)``, ``S = integral(r dl)`` and ``A_planar`` the plain
cell area. The distinction matters: the pressure source is integrated over the
planar area, not the r-weighted volume.

Note that the symmetry axis needs no special flux treatment. Faces lying on
r = 0 have ``S = 0``, so their flux contribution vanishes identically.
"""
import time

import numpy as np

from . import bc
from .reconstruct import get_limiter, muscl
from .riemann import get_flux
from .viscous import ViscousTerms

TINY = 1e-30

#: Low-storage multistage coefficients (Jameson). The 5-stage set is tuned for
#: steady-state convergence with local time stepping; rk3 is SSP and is what
#: you want for genuinely unsteady runs.
STAGE_COEFFS = {
    "rk1": (1.0,),
    "rk3": (0.1918, 0.4929, 1.0),
    "rk5": (0.0695, 0.1602, 0.2898, 0.5060, 1.0),
}


class NozzleSolver:
    """Steady-state (or time-accurate) axisymmetric nozzle solver.

    Parameters
    ----------
    grid : Grid
    gas : PerfectGas
    bcs : BoundaryConditions
    flux : {'roe', 'hllc', 'ausm'}
    limiter : {'vanalbada', 'minmod', 'vanleer', 'mc', 'none'}
        ``none`` gives a first-order scheme, useful for hard restarts.
    viscous : bool
        Set False for an Euler run (useful for verifying against the
        quasi-1-D isentropic solution).
    cfl : float
        Target CFL number. Ramped from ``cfl_start`` over ``cfl_ramp`` steps.
    visc_dt_factor : float
        Weight on the viscous spectral radius in the local time step.
    smooth_eps : float
        Implicit residual smoothing coefficient; 0 disables it. Around 0.5
        with two Jacobi sweeps permits roughly twice the CFL.
    """

    def __init__(self, grid, gas, bcs, flux="roe", limiter="vanalbada",
                 viscous=True, cfl=1.8, cfl_start=0.2, cfl_ramp=500,
                 entropy_fix=0.1, integrator="rk5", local_dt=True,
                 visc_dt_factor=2.0, smooth_eps=0.0, smooth_sweeps=2, ng=2):
        self.grid = grid
        self.gas = gas
        self.bcs = bcs
        self.ng = ng
        self.flux_fn = get_flux(flux)
        self.flux_name = flux
        self.limiter_fn = get_limiter(limiter)
        self.limiter_name = limiter
        self.viscous = bool(viscous)
        self.cfl = float(cfl)
        self.cfl_start = float(cfl_start)
        self.cfl_ramp = int(cfl_ramp)
        self.entropy_fix = float(entropy_fix)
        self.integrator = integrator
        self.local_dt = bool(local_dt)
        self.visc_dt_factor = float(visc_dt_factor)
        self.smooth_eps = float(smooth_eps)
        self.smooth_sweeps = int(smooth_sweeps)

        if integrator not in STAGE_COEFFS:
            raise ValueError(f"integrator must be one of {sorted(STAGE_COEFFS)}")

        self.vt = ViscousTerms(grid, gas, bcs, ng) if self.viscous else None

        ni, nj = grid.ni, grid.nj
        self.W = np.zeros((4, ni + 2 * ng, nj + 2 * ng))
        self.U = np.zeros((4, ni, nj))
        self.iter = 0
        self.history = []
        self._res_ref = None
        self.n_clip = 0

    # -- variable conversions --------------------------------------------
    def prim_to_cons(self, Wc):
        rho, u, v, p = Wc
        E = p / ((self.gas.gamma - 1.0) * rho) + 0.5 * (u * u + v * v)
        return np.stack([rho, rho * u, rho * v, rho * E])

    def cons_to_prim(self, U):
        rho = np.maximum(U[0], 1e-12)
        u = U[1] / rho
        v = U[2] / rho
        e = U[3] / rho - 0.5 * (u * u + v * v)
        p = (self.gas.gamma - 1.0) * rho * e
        bad = (p <= 0.0) | (U[0] <= 0.0) | ~np.isfinite(p)
        if bad.any():
            self.n_clip += int(bad.sum())
            p = np.where(bad, 1e-6 * max(self.bcs.p_amb, 1.0), p)
        return np.stack([rho, u, v, p])

    # -- initialisation ---------------------------------------------------
    def initialize(self, W_cells):
        """Seed the interior from a primitive field of shape (4, ni, nj)."""
        self.U = self.prim_to_cons(W_cells)
        self.iter = 0
        self.history = []
        self._res_ref = None
        return self

    def _fill(self, U):
        ng, ni, nj = self.ng, self.grid.ni, self.grid.nj
        self.W[:, ng:ng + ni, ng:ng + nj] = self.cons_to_prim(U)
        bc.apply(self.W, self.grid, self.gas, self.bcs, ng, viscous=self.viscous)
        return self.W

    # -- residual ---------------------------------------------------------
    def residual(self, U):
        """Net flux out of every cell, R, such that V dU/dt = -R."""
        return self.residual_from_W(self._fill(U))

    def residual_from_W(self, W):
        """Residual from an already ghost-filled primitive array.

        Split out from :meth:`residual` so that the discrete operator can be
        exercised on a prescribed state without the boundary conditions
        rewriting it -- see the free-stream preservation test.
        """
        grid, ng = self.grid, self.ng
        ni, nj = grid.ni, grid.nj
        Fi, Fj, tau_tt = self.face_fluxes(W)

        R = (Fi[:, 1:, :] * grid.S_i[1:, :] - Fi[:, :-1, :] * grid.S_i[:-1, :]
             + Fj[:, :, 1:] * grid.S_j[:, 1:] - Fj[:, :, :-1] * grid.S_j[:, :-1])

        p_cell = W[3, ng:ng + ni, ng:ng + nj]
        R[2] -= (p_cell - tau_tt) * grid.A_planar
        return R

    def face_fluxes(self, W):
        """Total flux (inviscid minus viscous) on both face families.

        These are the fluxes the update actually uses, so integrating them is
        the only self-consistent way to measure mass flow and thrust. Summing
        cell-centred states over a station instead introduces a quadrature
        error of its own, which on this geometry shows up as a spurious few
        per cent of apparent mass-flow variation along the nozzle.
        """
        grid, gas, ng = self.grid, self.gas, self.ng
        ni, nj = grid.ni, grid.nj

        WL, WR = muscl(W, 1, ng, ni, self.limiter_fn)
        WL = WL[:, :, ng:ng + nj]
        WR = WR[:, :, ng:ng + nj]
        Fi = self.flux_fn(WL, WR, grid.nx_i, grid.nr_i, gas,
                          entropy_fix=self.entropy_fix)

        WLj, WRj = muscl(W, 2, ng, nj, self.limiter_fn)
        WLj = WLj[:, ng:ng + ni, :]
        WRj = WRj[:, ng:ng + ni, :]
        Fj = self.flux_fn(WLj, WRj, grid.nx_j, grid.nr_j, gas,
                          entropy_fix=self.entropy_fix)

        # Solid wall: impose the pressure flux directly rather than relying on
        # the Riemann solver to reproduce it from mirrored states.
        p_w = WLj[3, :, -1]
        Fj[0, :, -1] = 0.0
        Fj[1, :, -1] = p_w * grid.nx_j[:, -1]
        Fj[2, :, -1] = p_w * grid.nr_j[:, -1]
        Fj[3, :, -1] = 0.0

        if self.viscous:
            Fvi, Fvj, tau_tt = self.vt.compute(W)
            Fi = Fi - Fvi
            Fj = Fj - Fvj
        else:
            tau_tt = 0.0
        return Fi, Fj, tau_tt

    # -- time step --------------------------------------------------------
    def local_timestep(self, U, cfl):
        grid, gas = self.grid, self.gas
        rho, u, v, p = self.cons_to_prim(U)
        a = np.sqrt(gas.gamma * p / rho)

        SI2 = grid.Sx_I ** 2 + grid.Sr_I ** 2
        SJ2 = grid.Sx_J ** 2 + grid.Sr_J ** 2
        lc_I = np.abs(u * grid.Sx_I + v * grid.Sr_I) + a * np.sqrt(SI2)
        lc_J = np.abs(u * grid.Sx_J + v * grid.Sr_J) + a * np.sqrt(SJ2)
        denom = lc_I + lc_J

        if self.viscous:
            T = p / (rho * gas.R)
            mu = gas.viscosity(T)
            c = np.maximum(4.0 / (3.0 * rho), gas.gamma / rho) * mu / gas.Pr
            lv_I = c * SI2 / grid.V
            lv_J = c * SJ2 / grid.V
            denom = denom + self.visc_dt_factor * (lv_I + lv_J)

        dt = cfl * grid.V / np.maximum(denom, TINY)
        if not self.local_dt:
            dt = np.full_like(dt, dt.min())
        return dt

    # -- residual smoothing -----------------------------------------------
    def _smooth(self, R):
        eps = self.smooth_eps
        if eps <= 0.0:
            return R
        Rs = R.copy()
        for _ in range(self.smooth_sweeps):
            acc = np.zeros_like(Rs)
            acc[:, :-1, :] += Rs[:, 1:, :]
            acc[:, 1:, :] += Rs[:, :-1, :]
            acc[:, :, :-1] += Rs[:, :, 1:]
            acc[:, :, 1:] += Rs[:, :, :-1]
            cnt = np.zeros(R.shape[1:])
            cnt[:-1, :] += 1
            cnt[1:, :] += 1
            cnt[:, :-1] += 1
            cnt[:, 1:] += 1
            Rs = (R + eps * acc) / (1.0 + eps * cnt)
        return Rs

    # -- stepping ---------------------------------------------------------
    def step(self):
        cfl = self.cfl
        if self.cfl_ramp > 0 and self.iter < self.cfl_ramp:
            f = self.iter / self.cfl_ramp
            cfl = self.cfl_start + f * (self.cfl - self.cfl_start)

        dt = self.local_timestep(self.U, cfl)
        fac = dt / self.grid.V

        U0 = self.U
        U = U0
        R = None
        for alpha in STAGE_COEFFS[self.integrator]:
            R = self.residual(U)
            U = U0 - alpha * fac * self._smooth(R)
        self.U = U
        self.iter += 1
        return R, dt

    def residual_norms(self, R):
        """RMS of dU/dt over the field, per equation."""
        d = R / self.grid.V
        n = d[0].size
        return np.sqrt(np.sum(d * d, axis=(1, 2)) / n)

    def run(self, max_iter=20000, tol=1e-6, print_every=250, callback=None):
        t0 = time.time()
        if not self.history:
            print(f"{'iter':>7} {'res(rho)':>12} {'drop':>10} "
                  f"{'dt_min':>11} {'CFL':>6} {'wall[s]':>9}")
        for _ in range(int(max_iter)):
            R, dt = self.step()
            norms = self.residual_norms(R)
            r = float(norms[0])
            if self._res_ref is None or not np.isfinite(self._res_ref) or self._res_ref == 0.0:
                self._res_ref = r if r > 0 else 1.0
            drop = r / self._res_ref
            self.history.append((self.iter, r, drop))

            if not np.isfinite(r):
                print(f"  diverged at iteration {self.iter}")
                return self

            if print_every and (self.iter % print_every == 0 or self.iter == 1):
                cfl = self.cfl if self.iter >= self.cfl_ramp else \
                    self.cfl_start + (self.iter / max(self.cfl_ramp, 1)) * (self.cfl - self.cfl_start)
                print(f"{self.iter:>7d} {r:>12.4e} {drop:>10.3e} "
                      f"{dt.min():>11.3e} {cfl:>6.2f} {time.time() - t0:>9.1f}")
            if callback is not None:
                callback(self)
            if drop < tol:
                print(f"  converged: residual dropped {drop:.2e} "
                      f"in {self.iter} iterations ({time.time() - t0:.1f} s)")
                break
        else:
            print(f"  stopped at {self.iter} iterations, residual drop "
                  f"{self.history[-1][2]:.2e} ({time.time() - t0:.1f} s)")
        if self.n_clip:
            print(f"  note: pressure floor applied {self.n_clip} cell-visits")
        return self

    # -- output -----------------------------------------------------------
    def primitives(self):
        """Interior primitive fields as (rho, u, v, p), each (ni, nj)."""
        return self.cons_to_prim(self.U)

    def fields(self):
        rho, u, v, p = self.primitives()
        gas = self.gas
        T = p / (rho * gas.R)
        a = np.sqrt(gas.gamma * p / rho)
        q = np.hypot(u, v)
        return dict(rho=rho, u=u, v=v, p=p, T=T, a=a, M=q / a, speed=q,
                    x=self.grid.xc, r=self.grid.rc)

    def save(self, path):
        np.savez_compressed(
            path, U=self.U, x_n=self.grid.x_n, r_n=self.grid.r_n,
            iter=self.iter, history=np.array(self.history))

    def load(self, path):
        d = np.load(path)
        if d["U"].shape != self.U.shape:
            raise ValueError("checkpoint grid does not match this solver")
        self.U = d["U"]
        self.iter = int(d["iter"])
        self.history = [tuple(row) for row in d["history"]]
        self._res_ref = self.history[0][1] if self.history else None
        return self
