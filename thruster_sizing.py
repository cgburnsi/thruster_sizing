# Thruster Sizing Code
import numpy as np
import metronos as mt


def _central_diff(f, x, h_rel=1e-6):
    h = abs(x) * h_rel + 1e-12
    return (f(x + h) - f(x - h)) / (2 * h)


class MachSolver:
    def __init__(self, tol=1e-4, max_iter=200):
        self.tol = tol
        self.max_iter = max_iter

    def Mach_f(self, k, expan, mach):
        term1 = 2 / (k + 1)
        term2 = (k - 1) / 2 * mach**2
        term3 = (term1 * (1 + term2))**((k + 1) / (2 * (k - 1)))
        return expan - (1 / mach) * term3

    def solve(self, guess=5.0, k=1.4, expan=4.0):
        mach = guess
        for _ in range(self.max_iter):
            f        = self.Mach_f(k, expan, mach)
            df       = _central_diff(lambda m: self.Mach_f(k, expan, m), mach)
            mach_new = mach - f / df

            if abs(f) < self.tol or abs(mach_new - mach) < self.tol:
                return mach_new

            mach = mach_new

        raise ValueError("Mach solution did not converge within the maximum number of iterations.")


class PcSolver:
    def __init__(self, tol=1e-4, max_iter=200):
        self.tol = tol
        self.max_iter = max_iter

    def Pc_f(self, k, expan, ambient_press, F_design, throat_area, exit_Mach, guess):
        a1 = F_design / (throat_area * k)
        a2 = 2 / (k - 1)
        a3 = 2 / (k + 1)
        a4 = (k + 1) / (k - 1)
        a5 = (k - 1) / k
        a6 = (k - 1) / 2
        a7 = k / (k - 1)

        a8 = 1 / (1 + a6 * exit_Mach**2)**a7
        a9 = np.sqrt(a2 * a3**a4 * (1 - a8**a5))

        return guess - (a1 * 1 / (a9 + expan * (a8 - (ambient_press / guess))))

    def solve(self, k, expan, ambient_press, F_design, throat_area, exit_Mach, guess=6e5):
        Pc = guess
        for _ in range(self.max_iter):
            f      = self.Pc_f(k, expan, ambient_press, F_design, throat_area, exit_Mach, Pc)
            df     = _central_diff(
                         lambda p: self.Pc_f(k, expan, ambient_press, F_design, throat_area, exit_Mach, p),
                         Pc)
            Pc_new = Pc - f / df

            if abs(f) < self.tol or abs(Pc_new - Pc) < self.tol:
                return Pc_new

            Pc = Pc_new

        raise ValueError("Pc solution did not converge within the maximum number of iterations.")


if __name__ == '__main__':

    # ── Inputs ───────────────────────────────────────────────────────────────
    d_c      = mt.Quantity(8.3,    'mm')   # Chamber diameter
    d_t      = mt.Quantity(0.29,   'mm')   # Throat diameter
    expan    = 100                          # Area expansion ratio (A_e/A_t)
    d_e      = d_t * np.sqrt(expan)        # Exit diameter
    F_set    = mt.Quantity(100,    'mN')   # Design thrust
    Pc_guess = mt.Quantity(100,    'psi')  # Chamber pressure initial guess
    P_a      = mt.Quantity(5,      'Torr') # Ambient pressure

    # Combustion gas properties (LOX/LH2, O/F = 5.0)
    eta         = 0.95
    k           = 1.26
    cstar_ideal = mt.Quantity(2350.0, 'm/s')
    MW          = 11.8
    T0          = mt.Quantity(3250.0, 'K')
    cstar       = cstar_ideal * eta

    # ── Geometry ─────────────────────────────────────────────────────────────
    A_c         = np.pi * d_c.to('m').value**2 / 4
    A_t         = np.pi * d_t.to('m').value**2 / 4
    A_e         = np.pi * d_e.to('m').value**2 / 4
    expan_ratio = A_e / A_t

    # ── Solve ─────────────────────────────────────────────────────────────────
    exit_Mach     = MachSolver().solve(k=k, expan=expan_ratio, guess=5.0)
    chamber_press = PcSolver().solve(k=k, expan=expan_ratio,
                                     ambient_press=P_a.to('Pa').value,
                                     F_design=F_set.to('N').value,
                                     throat_area=A_t,
                                     exit_Mach=exit_Mach,
                                     guess=Pc_guess.to('Pa').value)

    exit_press = chamber_press / (1 + ((k - 1) / 2) * exit_Mach**2)**(k / (k - 1))
    m_dot      = chamber_press * A_t / cstar.to('m/s').value
    Isp        = F_set.to('N').value / 9.807 / m_dot

    # ── Display ───────────────────────────────────────────────────────────────
    print("\nThruster Performance and Design Summary")
    print("-" * 80)
    print(f"{'Parameter':<25} {'SI Value':>15} {'Unit':<10} {'Eng Value':>15} {'Unit':<10}")
    print("-" * 80)

    print(f"{'Chamber Diameter, d_c':<25} {d_c.to('mm').value:>15.4f} {'mm':<10} {d_c.to('in').value:>15.4f} {'in':<10}")
    print(f"{'Throat Diameter, d_t':<25} {d_t.to('mm').value:>15.4f} {'mm':<10} {d_t.to('in').value:>15.4f} {'in':<10}")
    print(f"{'Exit Diameter, d_e':<25} {d_e.to('mm').value:>15.4f} {'mm':<10} {d_e.to('in').value:>15.4f} {'in':<10}")
    print(f"{'Chamber Area, A_c':<25} {mt.Quantity(A_c, 'm^2').to('mm^2').value:>15.4f} {'mm^2':<10} {mt.Quantity(A_c, 'm^2').to('in^2').value:>15.4f} {'in^2':<10}")
    print(f"{'Throat Area, A_t':<25} {mt.Quantity(A_t, 'm^2').to('mm^2').value:>15.4f} {'mm^2':<10} {mt.Quantity(A_t, 'm^2').to('in^2').value:>15.4f} {'in^2':<10}")
    print(f"{'Exit Area, A_e':<25} {mt.Quantity(A_e, 'm^2').to('mm^2').value:>15.4f} {'mm^2':<10} {mt.Quantity(A_e, 'm^2').to('in^2').value:>15.4f} {'in^2':<10}")
    print(f"{'Expansion Ratio':<25} {expan_ratio:>15.4f} {'-':<10} {'-':>15} {'-':<10}")
    print(f"{'Chamber Pressure, P_c':<25} {mt.Quantity(chamber_press, 'Pa').to('bar').value:>15.2f} {'bara':<10} {mt.Quantity(chamber_press, 'Pa').to('psi').value:>15.2f} {'psia':<10}")
    print(f"{'Exit Pressure, P_e':<25} {mt.Quantity(exit_press, 'Pa').to('bar').value:>15.4f} {'bara':<10} {mt.Quantity(exit_press, 'Pa').to('Torr').value:>15.2f} {'Torr':<10}")
    print(f"{'Ambient Pressure, P_a':<25} {P_a.to('bar').value:>15.4f} {'bara':<10} {P_a.to('Torr').value:>15.2f} {'Torr':<10}")
    print(f"{'Chamber Temperature, T_c':<25} {T0.to('degC').value:>15.2f} {'degC':<10} {T0.to('degF').value:>15.2f} {'degF':<10}")
    print(f"{'Exit Mach Number, M_e':<25} {exit_Mach:>15.4f} {'-':<10} {'-':>15} {'-':<10}")
    print(f"{'C* Efficiency, eta':<25} {eta * 100:>15.2f} {'%':<10} {'-':>15} {'-':<10}")
    print(f"{'Thrust, F':<25} {F_set.to('mN').value:>15.4f} {'mN':<10} {F_set.to('lbf').value:>15.4f} {'lbf':<10}")
    print(f"{'Mass Flowrate, mdot':<25} {mt.Quantity(m_dot, 'kg/s').to('g/s').value:>15.4f} {'g/s':<10} {mt.Quantity(m_dot, 'kg/s').to('lb/s').value:>15.4f} {'lb/s':<10}")
    print(f"{'Specific Impulse, Isp':<25} {Isp:>15.4f} {'sec':<10}  {'-':>15} {'-':<10}")
    print("-" * 80)

    # ── Catalyst bed sizing ───────────────────────────────────────────────────
    A_target   = 0.75
    d_foam     = mt.Quantity(8.3, 'mm').to('m').value
    L_foam     = mt.Quantity(20,  'mm').to('m').value
    A_foam     = mt.Quantity(np.pi * 8.3**2 / 4, 'mm^2').to('m^2').value
    V_foam     = A_foam * L_foam
    V_foam_cm3 = mt.Quantity(V_foam, 'm^3').to('cm^3').value

    PPI        = 60
    S_v        = 1300 * (PPI / 20)**1.3   # [m^2/m^3]
    A_total    = S_v * V_foam
    S_v_target = A_target / V_foam

    mdot       = mt.Quantity(0.4713, 'g/s').to('kg/s').value
    loading    = A_foam / mdot

    # ── Chamber wall thickness ────────────────────────────────────────────────
    yield_strength        = mt.Quantity(100, 'MPa').to('Pa').value
    pressure_maximum      = mt.Quantity(275, 'psi').to('Pa').value
    chamber_mean_diameter = 9.0
    wall_t = pressure_maximum * chamber_mean_diameter / (2 * yield_strength)
