# pfr_level0a.py
# Level 0a: Species-only PFR (isothermal, no dispersion), concentration form.
# Governing eq.: v * dc/dz = R(c) = nu @ r(c), with inlet c(0) = c_in.

import numpy as np

def rk4_step(f, z, y, dz):
    k1 = f(z, y)
    k2 = f(z + 0.5*dz, y + 0.5*dz*k1)
    k3 = f(z + 0.5*dz, y + 0.5*dz*k2)
    k4 = f(z + dz,     y + dz*k3)
    return y + (dz/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def integrate_pfr(v, z_end, nz, c_in, nu, kinetics_fn, clip_nonneg=True):
    """
    v           : superficial velocity [m/s] (constant)
    z_end       : reactor length [m]
    nz          : number of axial nodes (uniform grid)
    c_in        : 1D array [ns] of inlet concentrations [mol/m^3]
    nu          : stoichiometry matrix shape [ns, nr] (prod +, react -)
    kinetics_fn : function r = kinetics_fn(c)
                  returns reaction-rate vector r shape [nr] in [mol_rxn/m^3/s]
    clip_nonneg : if True, enforce c>=0 after each step

    returns z (centers), C profile array shape [nz, ns]
    """
    c_in = np.asarray(c_in, dtype=float)
    nu   = np.asarray(nu, dtype=float)

    ns = c_in.size
    assert nu.shape[0] == ns, "nu must have ns rows"

    # Grid
    z_faces = np.linspace(0.0, z_end, nz + 1)
    dz      = z_faces[1] - z_faces[0]
    zc      = 0.5*(z_faces[:-1] + z_faces[1:])  # cell centers

    C = np.zeros((nz, ns), dtype=float)
    c = c_in.copy()

    # RHS as function of z and c: dc/dz = (nu @ r(c)) / v
    def rhs(_, c_local):
        r = np.asarray(kinetics_fn(c_local), dtype=float)  # [nr]
        R = nu.dot(r)                                      # [ns]
        return R / max(v, 1e-30)

    # March cell-by-cell with RK4 using face-to-face steps
    for k in range(nz):
        c = rk4_step(rhs, z_faces[k], c, dz)
        if clip_nonneg:
            c = np.maximum(c, 0.0)
        C[k, :] = c

    return zc, C

# ---------------------------
# Example: A -> B, first-order in A
# ---------------------------
def example_first_order():
    # Species order: [A, B]
    k = 0.8      # 1/s
    v = 0.2      # m/s
    L = 1.0      # m
    nz = 200

    cA_in = 5.0  # mol/m^3
    cB_in = 0.0

    # Stoichiometry for A -> B: nu = [[-1], [+1]] (ns=2, nr=1)
    nu = np.array([[-1.0],
                   [ +1.0]])

    def kinetics(c):
        cA = c[0]
        r1 = k * max(cA, 0.0)   # mol_rxn/m^3/s
        return np.array([r1])

    z, C = integrate_pfr(v=v, z_end=L, nz=nz, c_in=[cA_in, cB_in], nu=nu, kinetics_fn=kinetics)
    cA = C[:, 0]; cB = C[:, 1]

    # Sanity: analytical solution for first-order A -> B
    cA_analytical = cA_in * np.exp(-(k/v) * z)
    max_rel_err = np.max(np.abs(cA - cA_analytical) / np.maximum(cA_analytical, 1e-12))
    X = 1.0 - cA[-1] / cA_in

    print(f"Outlet: c_A = {cA[-1]:.6f} mol/m^3, c_B = {cB[-1]:.6f} mol/m^3")
    print(f"Conversion of A = {100*X:.3f}%")
    print(f"Max relative error vs analytical = {100*max_rel_err:.4f}%")

    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(z, cA, label='c_A (num)')
        plt.plot(z, cB, label='c_B (num)')
        plt.plot(z, cA_analytical, '--', label='c_A (analytical)')
        plt.xlabel('z [m]'); plt.ylabel('Concentration [mol/m^3]')
        plt.title('Level 0a: A -> B, first-order')
        plt.legend(); plt.tight_layout(); plt.show()
    except Exception:
        pass

# ---------------------------
# Extension hooks (keep for later levels)
# ---------------------------
# - To switch the independent variable to reactor volume V = A z:
#     dc/dV = R(c) / F_v  ; you’d integrate with dV instead of dz.
# - To add more reactions: supply a larger nu [ns, nr] and return r(c) with length nr.
# - To add temperature later: pass T through kinetics_fn or close with enthalpy.

if __name__ == '__main__':
    example_first_order()


