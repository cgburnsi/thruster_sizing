"""Viscous fluxes and the axisymmetric hoop-stress source term.

Cell-centred gradients come from a Green-Gauss loop over the planar cell
faces. Face gradients are then formed by averaging the two adjoining cell
gradients and applying the standard directional correction, which removes the
odd-even decoupling that a plain average suffers on stretched grids.

The wall is handled exactly rather than by averaging: because ``u = v = 0``
everywhere along the wall, the velocity gradient there is purely wall-normal,
so it can be built directly from the first cell centre and its wall distance.
"""
import numpy as np

TINY = 1e-30


def _green_gauss(phi_if, phi_jf, grid):
    """Planar cell gradient from face values. Returns (d/dx, d/dr)."""
    gx = (phi_if[1:, :] * grid.nx_i[1:, :] * grid.L_i[1:, :]
          - phi_if[:-1, :] * grid.nx_i[:-1, :] * grid.L_i[:-1, :]
          + phi_jf[:, 1:] * grid.nx_j[:, 1:] * grid.L_j[:, 1:]
          - phi_jf[:, :-1] * grid.nx_j[:, :-1] * grid.L_j[:, :-1])
    gr = (phi_if[1:, :] * grid.nr_i[1:, :] * grid.L_i[1:, :]
          - phi_if[:-1, :] * grid.nr_i[:-1, :] * grid.L_i[:-1, :]
          + phi_jf[:, 1:] * grid.nr_j[:, 1:] * grid.L_j[:, 1:]
          - phi_jf[:, :-1] * grid.nr_j[:, :-1] * grid.L_j[:, :-1])
    return gx / grid.A_planar, gr / grid.A_planar


def _face_avg(phi, ng, ni, nj):
    """Face-averaged values on both face families, from a ghosted cell array."""
    js = slice(ng, ng + nj)
    isl = slice(ng, ng + ni)
    f_i = 0.5 * (phi[ng - 1:ng + ni, js] + phi[ng:ng + ni + 1, js])
    f_j = 0.5 * (phi[isl, ng - 1:ng + nj] + phi[isl, ng:ng + nj + 1])
    return f_i, f_j


def _corrected(gL_x, gL_r, gR_x, gR_r, phiL, phiR, dx, dr, dist):
    ex, er = dx / dist, dr / dist
    ax = 0.5 * (gL_x + gR_x)
    ar = 0.5 * (gL_r + gR_r)
    corr = (phiR - phiL) / dist - (ax * ex + ar * er)
    return ax + corr * ex, ar + corr * er


class ViscousTerms:
    """Assembles viscous face fluxes and the hoop-stress source."""

    def __init__(self, grid, gas, bcs, ng=2):
        self.grid = grid
        self.gas = gas
        self.bcs = bcs
        self.ng = ng
        self.dxi, self.dri, self.disti = grid.cell_centre_delta_i()
        self.dxj, self.drj, self.distj = grid.cell_centre_delta_j()
        # Normal distance from the last cell centre to the wall face midpoint
        self.d_wall = np.abs((grid.xf_j[:, -1] - grid.xc[:, -1]) * grid.nx_j[:, -1]
                             + (grid.rf_j[:, -1] - grid.rc[:, -1]) * grid.nr_j[:, -1])
        self.d_wall = np.maximum(self.d_wall, TINY)

    def compute(self, W):
        grid, gas, ng = self.grid, self.gas, self.ng
        ni, nj = grid.ni, grid.nj
        js = slice(ng, ng + nj)
        isl = slice(ng, ng + ni)

        rho, u, v, p = W
        T = p / (rho * gas.R)

        u_if, u_jf = _face_avg(u, ng, ni, nj)
        v_if, v_jf = _face_avg(v, ng, ni, nj)
        T_if, T_jf = _face_avg(T, ng, ni, nj)

        # Cell-centred gradients
        gux, gur = _green_gauss(u_if, u_jf, grid)
        gvx, gvr = _green_gauss(v_if, v_jf, grid)
        gTx, gTr = _green_gauss(T_if, T_jf, grid)

        uc = u[isl, js]
        vc = v[isl, js]
        Tc = T[isl, js]

        # ---- i-face gradients -------------------------------------------
        Gux_i = np.empty((ni + 1, nj))
        Gur_i = np.empty((ni + 1, nj))
        Gvx_i = np.empty((ni + 1, nj))
        Gvr_i = np.empty((ni + 1, nj))
        GTx_i = np.empty((ni + 1, nj))
        GTr_i = np.empty((ni + 1, nj))

        for (Gx, Gr, gx, gr, c) in (
                (Gux_i, Gur_i, gux, gur, uc),
                (Gvx_i, Gvr_i, gvx, gvr, vc),
                (GTx_i, GTr_i, gTx, gTr, Tc)):
            fx, fr = _corrected(gx[:-1], gr[:-1], gx[1:], gr[1:],
                                c[:-1], c[1:], self.dxi, self.dri, self.disti)
            Gx[1:-1] = fx
            Gr[1:-1] = fr
            Gx[0], Gr[0] = gx[0], gr[0]        # inlet: use the adjacent cell
            Gx[-1], Gr[-1] = gx[-1], gr[-1]    # outlet: likewise

        # ---- j-face gradients -------------------------------------------
        Gux_j = np.empty((ni, nj + 1))
        Gur_j = np.empty((ni, nj + 1))
        Gvx_j = np.empty((ni, nj + 1))
        Gvr_j = np.empty((ni, nj + 1))
        GTx_j = np.empty((ni, nj + 1))
        GTr_j = np.empty((ni, nj + 1))

        for (Gx, Gr, gx, gr, c) in (
                (Gux_j, Gur_j, gux, gur, uc),
                (Gvx_j, Gvr_j, gvx, gvr, vc),
                (GTx_j, GTr_j, gTx, gTr, Tc)):
            fx, fr = _corrected(gx[:, :-1], gr[:, :-1], gx[:, 1:], gr[:, 1:],
                                c[:, :-1], c[:, 1:], self.dxj, self.drj, self.distj)
            Gx[:, 1:-1] = fx
            Gr[:, 1:-1] = fr
            Gx[:, 0], Gr[:, 0] = gx[:, 0], gr[:, 0]   # axis face carries zero area

        # Wall face: gradient is exactly normal because u = v = 0 along the wall
        nxw, nrw = grid.nx_j[:, -1], grid.nr_j[:, -1]
        dudn = -uc[:, -1] / self.d_wall
        dvdn = -vc[:, -1] / self.d_wall
        Gux_j[:, -1], Gur_j[:, -1] = dudn * nxw, dudn * nrw
        Gvx_j[:, -1], Gvr_j[:, -1] = dvdn * nxw, dvdn * nrw
        if self.bcs.wall_is_adiabatic:
            dTdn = np.zeros_like(dudn)
        else:
            dTdn = (self.bcs.T_wall - Tc[:, -1]) / self.d_wall
        GTx_j[:, -1], GTr_j[:, -1] = dTdn * nxw, dTdn * nrw

        # ---- assemble the fluxes ----------------------------------------
        Fv_i = self._flux(u_if, v_if, T_if, grid.rf_i,
                          Gux_i, Gur_i, Gvx_i, Gvr_i, GTx_i, GTr_i,
                          grid.nx_i, grid.nr_i)
        Fv_j = self._flux(u_jf, v_jf, T_jf, grid.rf_j,
                          Gux_j, Gur_j, Gvx_j, Gvr_j, GTx_j, GTr_j,
                          grid.nx_j, grid.nr_j)

        # ---- hoop stress (cell-centred), used by the axisymmetric source --
        mu = gas.viscosity(Tc)
        divV = gux + gvr + vc / np.maximum(grid.rc, TINY)
        tau_tt = 2.0 * mu * vc / np.maximum(grid.rc, TINY) - (2.0 / 3.0) * mu * divV
        return Fv_i, Fv_j, tau_tt

    def _flux(self, uf, vf, Tf, rf, gux, gur, gvx, gvr, gTx, gTr, nx, nr):
        gas = self.gas
        mu = gas.viscosity(Tf)
        kap = gas.conductivity(Tf)

        # v/r degenerates on the axis, where symmetry gives v/r -> dv/dr
        vr = np.where(rf > TINY, vf / np.where(rf > TINY, rf, 1.0), gvr)
        divV = gux + gvr + vr
        lam = -(2.0 / 3.0) * mu

        txx = 2.0 * mu * gux + lam * divV
        trr = 2.0 * mu * gvr + lam * divV
        txr = mu * (gur + gvx)
        qx = -kap * gTx
        qr = -kap * gTr

        fx = txx * nx + txr * nr
        fr = txr * nx + trr * nr
        fe = ((uf * txx + vf * txr - qx) * nx + (uf * txr + vf * trr - qr) * nr)
        return np.stack([np.zeros_like(fx), fx, fr, fe])
