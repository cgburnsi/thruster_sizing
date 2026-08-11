# pfr_level0a_solver.py
# Minimal PFR solver class: species-only, isothermal, no dispersion (Level 0a).
# Governing eq.: v * dc/dz = nu @ r(c), with inlet c(0) = c_in.

import numpy as np

class PFRSolver:
    def __init__(self, v, L, nz, c_in, nu, kinetics_fn, clip_nonneg=True):
        """
        v           : superficial velocity [m/s] (constant)
        L           : reactor length [m]
        nz          : number of axial cells (uniform)
        c_in        : array-like [ns] inlet concentrations [mol/m^3]
        nu          : stoichiometry matrix [ns, nr] (prod +, react -)
        kinetics_fn : function r = kinetics_fn(c) -> [nr] [mol_rxn/m^3/s]
        clip_nonneg : if True, enforce c>=0 after each step
        """
        self.v = float(v)
        self.L = float(L)
        self.nz = int(nz)
        self.c_in = np.asarray(c_in, dtype=float).copy()
        self.nu = np.asarray(nu, dtype=float).copy()
        self.kinetics_fn = kinetics_fn
        self.clip_nonneg = bool(clip_nonneg)

        assert self.nu.shape[0] == self.c_in.size, "nu must have ns rows"
        self.ns = self.c_in.size
        self.nr = self.nu.shape[1]

        # Grid (face-to-face marching; store centers for output)
        self.z_faces = np.linspace(0.0, self.L, self.nz + 1)
        self.dz = self.z_faces[1] - self.z_faces[0]
        self.z = 0.5 * (self.z_faces[:-1] + self.z_faces[1:])

        # Storage for last run
        self.C = None

    def rhs(self, z, c):
        """dc/dz = (nu @ r(c)) / v"""
        r = np.asarray(self.kinetics_fn(c), dtype=float)
        if r.shape != (self.nr,):
            raise ValueError("kinetics_fn must return length nr")
        R = self.nu.dot(r)   # [ns]
        return R / max(self.v, 1e-30)

    def step_rk4(self, z, c, dz):
        k1 = self.rhs(z, c)
        k2 = self.rhs(z + 0.5*dz, c + 0.5*dz*k1)
        k3 = self.rhs(z + 0.5*dz, c + 0.5*dz*k2)
        k4 = self.rhs(z + dz,     c + dz*k3)
        return c + (dz/6.0)*(k1 + 2*k2 + 2*k3 + k4)

    def solve(self, method='rk4'):
        """Integrate from z=0 to L; returns (z_centers, C[nz, ns])."""
        c = self.c_in.copy()
        C = np.zeros((self.nz, self.ns), dtype=float)

        step = self.step_rk4 if method.lower() == 'rk4' else None
        if step is None:
            raise ValueError("Only 'rk4' is implemented in this minimal solver")

        for k in range(self.nz):
            c = step(self.z_faces[k], c, self.dz)
            if self.clip_nonneg:
                c = np.maximum(c, 0.0)
            C[k, :] = c

        self.C = C
        return self.z, C

    # ---- Convenience utilities ----
    def outlet(self):
        if self.C is None:
            raise RuntimeError("Run solve() first")
        return self.C[-1, :].copy()

    def conversions(self, which):
        """Return conversion of species index 'which' against its inlet."""
        if self.C is None:
            raise RuntimeError("Run solve() first")
        c0 = max(self.c_in[which], 1e-30)
        return 1.0 - self.C[-1, which] / c0

    def replace_kinetics(self, kinetics_fn):
        self.kinetics_fn = kinetics_fn

    # ---- Extension hooks you can add later ----
    # - allow v(z) as a callable
    # - dispersion term with Danckwerts BCs
    # - enthalpy energy equation and coupling
    # - volume-based marching (d/dV) instead of d/dz

# ---------------------------
# Example: A + B -> C (second order)
# ---------------------------
def _example():
    # Species order: [A, B, C]
    k2 = 0.10     # m^3/(mol*s)
    v  = 0.20     # m/s
    L  = 2.0      # m
    nz = 400

    c_in = [5.0, 3.0, 0.0]  # mol/m^3

    # Stoichiometry: A + B -> C  ==> nu = [-1, -1, +1]^T
    nu = np.array([[-1.0],
                   [-1.0],
                   [+1.0]])

    def kinetics(c):
        cA, cB, cC = c
        r1 = k2 * max(cA, 0.0) * max(cB, 0.0)  # mol_rxn/m^3/s
        return np.array([r1])

    solver = PFRSolver(v=v, L=L, nz=nz, c_in=c_in, nu=nu, kinetics_fn=kinetics)
    z, C = solver.solve()
    cA, cB, cC = C[:,0], C[:,1], C[:,2]

    cA_out, cB_out, cC_out = solver.outlet()
    XA = solver.conversions(which=0)
    XB = solver.conversions(which=1)

    print(f"Outlet: c_A={cA_out:.4f}, c_B={cB_out:.4f}, c_C={cC_out:.4f}  [mol/m^3]")
    print(f"Conversions: X_A={100*XA:.2f}%, X_B={100*XB:.2f}%")

    # Quick invariants for A+B->C
    inv1 = np.max(np.abs((cA - cB) - (c_in[0] - c_in[1])))
    inv2 = np.max(np.abs(cC - (c_in[2] + (c_in[0] - cA))))
    print(f"Invariants check -> max| (A-B) - (A0-B0) | = {inv1:.3e}, max| C - (C0 + extent) | = {inv2:.3e}")

    try:
        import matplotlib.pyplot as plt
        import matplotlib
        plt.figure(); plt.plot(z, cA, label='A'); plt.plot(z, cB, label='B'); plt.plot(z, cC, label='C')
        plt.xlabel('z [m]'); plt.ylabel('Concentration [mol/m^3]'); plt.title('PFR Level 0a: A + B → C'); plt.legend(); plt.tight_layout(); plt.show()
    except Exception:
        pass

if __name__ == '__main__':
    
    #_example()

     
    RES = XD, YD = 100, 100   
