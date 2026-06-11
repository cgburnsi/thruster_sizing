# Thruster Sizing Code
import convert as cv
import numpy as np

class MachSolver:
    def __init__(self, tol=1e-4, max_iter=200):
        self.tol = tol                           # [-] Convergence Tolerance
        self.max_iter = max_iter                 # [-] Maximum allowable iterations to reach convergence

    def Mach_f(self, k, expan, mach):
        term1   = 2 / (k + 1)
        term2   = (k - 1) / 2 * mach**2
        term3   = (term1 * (1 + term2))**((k + 1) / (2 * (k - 1)))
    
        return expan - (1 / mach) * term3
    
    def Mach_df(self, k, mach):
        a1 = k + 1
        a2 = k - 1
        a3 = a1 / a2
        a4 = 2 / a1
        a5 = a2 / 2
        a6 = a3 / 2
        
        term1 = (a4 * (a5 * mach**2 + 1))**a6 / mach**2
        term2 = a3 * a4 * a5 * (a4 * (a5 * mach**2 + 1))**(a6 - 1)
        
        return term1 - term2


    def solve(self, guess=5.0, k=1.4, expan=4.0):
        mach = guess  # Initialize with the user's guess
        for _ in range(self.max_iter):
            f           = self.Mach_f(k, expan, mach)   # Mach Number - Area Ratio Calculation
            df          = self.Mach_df(k, mach)         # Analytical Derivative of the Mach Number-Area Ratio    
            mach_new    = mach - f / df                 # Update mach with the new value
    
            # Convergence check: either function value or change in mach number
            if abs(f) < self.tol or abs(mach_new - mach) < self.tol:
                return mach_new
    
            mach        = mach_new                      # Update mach for the next iteration
        
        # If the loop ends without convergence
        raise ValueError("Solution did not converge within the maximum number of iterations.")


class PcSolver:
    def __init__(self, tol=1e-4, max_iter=200):
        self.tol = tol                 # Convergence tolerance
        self.max_iter = max_iter       # Maximum number of iterations

    def Pc_f(self, k, expan, ambient_press, F_design, throat_area, exit_Mach, guess):
        # Calculate the estimated chamber pressure using the isentropic 
        # relation, area ratio, and the coefficient of thrust        
        a1 = F_design / (throat_area * k)
        a2 = 2 / (k - 1)
        a3 = 2 / (k + 1)
        a4 = (k + 1) / (k - 1)
        a5 = (k - 1) / k
        a6 = (k - 1) / 2
        a7 = k / (k - 1)
        
        a8 = 1 / (1 + a6 * exit_Mach**2)**a7
        a9 = np.sqrt(a2 * a3**a4 * (1 - a8**a5))

        # The function output
        return guess - (a1 * 1 / (a9 + expan * (a8 - (ambient_press / guess))))

    def Pc_df(self, k, expan, ambient_press, F_design, throat_area, exit_Mach, guess):
        # Derivative of the function in 
        a1 = F_design * ambient_press * expan        
        df = a1/(throat_area*guess**2*k*(2**(1/2)*(-((2/(k + 1))**((k + 1)/(k - 1))* \
             ((1/((k/2 - 1/2)*exit_Mach**2 + 1)**(k/(k - 1)))**((k - 1)/k) - 1)) / \
             (k - 1))**(1/2) - expan*(ambient_press/guess - 1/((k/2 - 1/2)*exit_Mach**2 + 1)** \
             (k/(k - 1))))**2) + 1
        
        return df

    def solve(self, k, expan, ambient_press, F_design, throat_area, exit_Mach, guess=6e5):
        Pc = guess
        for _ in range(self.max_iter):
            f = self.Pc_f(k, expan, ambient_press, F_design, throat_area, exit_Mach, Pc)
            df = self.Pc_df(k, expan, ambient_press, F_design, throat_area, exit_Mach, Pc)
            Pc_new = Pc - f / df
            
            if abs(f) < self.tol or abs(Pc_new - Pc) < self.tol:
                return Pc_new
            
            Pc = Pc_new
        
        # If we exit the loop without convergence, raise an error
        raise ValueError("Solution did not converge within the maximum number of iterations.")


    

class Thruster:
    def __init__(self, chamber_dia, throat_dia, exit_dia, F_design, eta_cstar, Pc_guess, ambient_press):
        self.chamber_dia = chamber_dia
        self.throat_dia = throat_dia
        self.exit_dia = exit_dia
        self.F_design = F_design
        self.eta_cstar = eta_cstar
        self.Pc_guess = Pc_guess
        self.ambient_press  = ambient_press

        self.expan = (self.exit_dia / self.throat_dia)**2  # [-] Nozzle Expansion Ratio
        self.chamber_area = (np.pi * self.chamber_dia**2) / 4  # [m^2] Chamber Area
        self.throat_area = (np.pi * self.throat_dia**2) / 4    # [m^2] Throat Area
        self.exit_area = (np.pi * self.exit_dia**2) / 4        # [m^2] Exit Area

        # ASCENT Combustion Properties
        self.k = 1.2105  # [-] Combustion Chamber Ratio of Specific Heats
        self.cstar_ideal = 1396.0  # [m/s] Combustion Chamber ideal C*
        self.MW = 21.725  # [kg/kg-mol] Combustion Chamber Molecular Mass
        self.T0 = 2158.88  # [degK] Combustion Chamber Temperature
        self.cstar = self.cstar_ideal * self.eta_cstar  # [m/s] Estimated Chamber C*

        # Newton-Raphson Solver
        self.Mach   = MachSolver(tol=1e-4, max_iter=200)
        self.Pc     = PcSolver(tol=1e-4, max_iter=200)

        # Solve for exit Mach number
        self.exit_Mach      = self.Mach.solve(k=self.k, expan=self.expan, guess=5.0)
        self.chamber_press  = self.Pc.solve(k=self.k, expan=self.expan, ambient_press=self.ambient_press, F_design=self.F_design, 
                                            throat_area = self.throat_area, exit_Mach=self.exit_Mach, guess=6e5)
        
        self.exit_press = self.chamber_press /(1+((self.k-1)/2)*self.exit_Mach**2)**(self.k/(self.k-1))
        self.m_dot =  self.chamber_press * self.throat_area / self.cstar
        self.Isp =     self.F_design/9.807/self.m_dot
    
    def display_thruster_info_with_units(thruster):   
        print("\nThruster Performance and Design Summary")
        print("-" * 80)
        print(f"{'Parameter':<25} {'SI Value':>15} {'Unit':<10} {'Eng Value':>15} {'Unit':<10}")
        print("-" * 80)
    
        # Add each parameter with SI and English unit equivalents
        print(f"{'Chamber Diameter, d_c':<25} {cv.convert(thruster.chamber_dia, 'm', 'mm'):>15.4f} {'mm':<10} {cv.convert(thruster.chamber_dia, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Throat Diameter, d_t':<25} {cv.convert(thruster.throat_dia,'m', 'mm'):>15.4f} {'mm':<10} {cv.convert(thruster.throat_dia, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Exit Diameter, d_e':<25} {cv.convert(thruster.exit_dia, 'm', 'mm'):>15.4f} {'mm':<10} {cv.convert(thruster.exit_dia, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Chamber Area, A_c':<25} {cv.convert(thruster.chamber_area, 'm^2', 'mm^2'):>15.4f} {'mm^2':<10} {cv.convert(thruster.chamber_area, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Throat Area, A_t':<25} {cv.convert(thruster.throat_area, 'm^2', 'mm^2'):>15.4f} {'mm^2':<10} {cv.convert(thruster.throat_area, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Exit Area, A_e':<25} {cv.convert(thruster.exit_area, 'm^2', 'mm^2'):>15.4f} {'mm^2':<10} {cv.convert(thruster.exit_area, 'm', 'in'):>15.4f} {'in':<10}")
        print(f"{'Expansion Ratio':<25} {thruster.expan:>15.4f} {'-':<10} {'-':>15} {'-':<10}")
        print(f"{'Chamber Pressure, P_c':<25} {cv.convert(thruster.chamber_press, 'Pa', 'bar'):>15.2f} {'bara':<10} {cv.convert(thruster.chamber_press, 'Pa', 'psi'):>15.2f} {'psia':<10}")
        print(f"{'Exit Pressure, P_e':<25} {cv.convert(thruster.exit_press, 'Pa', 'bar'):>15.2f} {'bara':<10} {cv.convert(thruster.exit_press, 'Pa', 'Torr'):>15.2f} {'Torr':<10}")
        print(f"{'Ambient Pressure, P_a':<25} {cv.convert(thruster.ambient_press, 'Pa', 'bar'):>15.2f} {'bara':<10} {cv.convert(thruster.ambient_press, 'Pa', 'Torr'):>15.2f} {'Torr':<10}")
        print(f"{'Chamber Temperature, T_c':<25} {cv.convert(thruster.T0, 'K', 'degC'):>15.2f} {'degC':<10} {cv.convert(thruster.T0, 'K', 'degF'):>15.2f} {'degF':<10}")
        print(f"{'Exit Mach Number, M_e':<25} {thruster.exit_Mach:>15.4f} {'-':<10} {'-':>15} {'-':<10}")
        print(f"{'C* Efficiency, eta_cstar':<25} {thruster.eta_cstar * 100:>15.2f} {'%':<10} {'-':>15} {'-':<10}")
        print(f"{'Thrust, F':<25} {cv.convert(thruster.F_design, 'N', 'mN'):>15.4f} {'mN':<10} {cv.convert(thruster.F_design, 'N', 'lbf'):>15.4f} {'lbf':<10}")
        print(f"{'Mass Flowrate, mdot':<25} {cv.convert(thruster.m_dot, 'kg/s', 'g/s'):>15.4f} {'g/s':<10} {cv.convert(thruster.m_dot, 'kg/s', 'lb/s'):>15.4f} {'lb/s':<10}")
        print(f"{'Specific Impulse, Isp':<25} {thruster.Isp:>15.4f} {'sec':<10}  {'-':>15} {'-':<10}")
        print("-" * 80)

if __name__ == '__main__':
    
    # Geometry
    chamber_dia     = cv.convert(8.3, 'mm', 'm')    # [m] Chamber Diameter
    throat_dia      = cv.convert(.5, 'mm', 'm')    # [m] Throat Diameter
    exit_dia        = cv.convert(4.0, 'mm', 'm')    # [m] Exit Diameter

    # Performance
    F_design        = cv.convert(500, 'mN', 'N')    # [N] Thruster Thrust Design Point
    eta_cstar       = 0.92                           # [-] C* Efficiency
    Pc_guess        = cv.convert(100, 'psi', 'Pa')  # [Pa] Estimated Chamber Pressure
    
    # Operating Conditions
    ambient_press   = cv.convert(20, 'Torr', 'Pa')
    
    t = Thruster(chamber_dia, throat_dia, exit_dia, F_design, eta_cstar, Pc_guess, ambient_press)
    t.display_thruster_info_with_units()

    # Calculate catalyst Length
    A_target = 0.75              # [m^2] Catalytic Surface Area for 1N Thruster
    d_foam = cv.convert(8.3, 'mm', 'm')
    L_foam = cv.convert(20, 'mm', 'm')
    v_est = L_foam*(np.pi*d_foam**2)/4
    A_foam = cv.convert(np.pi *8.3**2/4, 'mm^2', 'm^2')   # [mm^2 Cross-Sectional Area of Foam
    V_foam = A_foam * cv.convert(20, 'mm', 'm')    
    V_foam_cm3 = cv.convert(V_foam, 'm^3', 'cm^3')
    
    PPI = 60
    
    S_v = 1300 * (PPI/20) ** (1.3)  # [m**2/m**3]
    A_total = S_v * V_foam
    A_total_cm2 = cv.convert(A_total, 'm^3', 'cm^3')

    S_v_target = A_target / V_foam


    mdot = cv.convert(0.4713, 'g/s', 'kg/s')
    
    loading = A_foam / mdot

# Estimate chamber wall thickness
yield_strength = cv.convert(100, 'MPa', 'Pa')
pressure_maximum = cv.convert(275, 'psi', 'Pa')
chamber_mean_diameter = 9.0

wall_t = pressure_maximum * chamber_mean_diameter / (2 * yield_strength)





