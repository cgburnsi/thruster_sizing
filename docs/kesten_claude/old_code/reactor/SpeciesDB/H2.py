# species/H2.py
import numpy as np
import matplotlib.pyplot as plt
from species.base_species import BaseSpecies

''' Fluid Properties - H2

    This class allows the user to use chemical species information useful for combustion 
    systems or fluid systems typically used in rocket propellant systems.  This class 
    (and others in this package) uses historical data verified to aerospace standards.
    A more modern method to determine fluid properties is to use CoolProp or REFPROP.
    Those codes will produce similar data and are generally more complete over wider 
    ranges of conditions, but the species list is limited.  
    
    Most historical datasets use English units.  This class uses SI units from simple 
    conversion of the English units.  The conversion is external to this class.  
    No provisions are provided for using English units.  If the user wishes to convert
    this class to English units, the values are provided in the comments.
    
    The code below uses the abbreviation Cp for specific heat.  The original Kesten
    paper and fortran code use the abbreviation Cf.  In the case of Kesten, this 
    was used to show the data was 'fitted' and not a purely thermodynamic source.  It
    is also typical for the time period when the code was written.
    
    If you investigate the plot of Cp vs. Temperature you will see a dip at 1000 K.  
    I've tried to find why this is the case, but i have not been successful.  This 
    may be an error in the Kesten data blocks, but I can't be sure.  The plot is
    not as smooth as I was expecting it to be.
    
    # SPECIFIC HEAT INFORMATION
        Source(s):
            Kesten, Arthur S., 'Analytical Study of Catalytic Reactors For Hydrazine
            Decomposition: One- and Two-Dimensional Steady-State Programs' 
            (Report No. NASA NOT AVAILABLE, Report No. UACRL G910461-30), 
            UACRL (United Aircraft Research Laboratories), August 1968
    
        # Specific Heat Data (BTU/lb-degR)
        T      = [540.0, 720.0, 900.0, 1080.0, 1260.0, 1440.0, 1620.0, 1800.0, 
                  1980.0, 2160.0, 2340.0, 2520.0, 2700.0, 2880.0, 3060.0]
        CFTBL1 = [3.4194, 3.4596, 3.4685, 3.4765, 3.4899, 3.5151, 3.5454, 3.5006,
                  3.6208, 3.6654, 3.7150, 3.7696, 3.8291, 3.8802, 3.9288]            

    # VISCOSITY INFORMATION
        Source(s):
            Kesten, Arthur S., 'Study of Catalytic Reactors for Hydrogen-Oxygen
            Ignition' (Report No. NASA CR-72567, Report No. UARL H910721), 
            UARL (United Aircraft Research Laboratories), July 1969
    
        # Viscosity Data (lb/ft-sec)
        T    = [180.0, 360.0, 720.0, 1080.0, 1440.0, 1800.0, 2160.0, 2520.0, 
                2880.0, 3240.0, 3600.0, 3960.0, 4320.0, 4680.0, 5040.0]
        VISC = [2.54e-6, 4.47e-6, 7.30e-6, 9.54e-6, 11.48e-6, 13.29e-6, 
                14.97e-6, 16.54e-6, 18.04e-6, 19.47e-6, 20.84e-6, 22.16e-6, 
                23.45e-6, 24.68e-6, 25.89e-6] 

    This class doesn't do any error catching.  For example, it will not throw an 
    error if you attempt to interpolate outside of the temperature ranges of 
    the data.  It will toss out np.NaN though, but you might get pretty far in 
    your calculations before it becomes a problem.
        
    TODO - Make a method to output a csv list or a pandas dataframe.  Yeah!
    TODO - Add boundary checks to the calculations/look-ups.  
           Maybe make it two functions.  One with and one without bound check
           so I can go faster if I ever had a problem with speed.
    
    Attributes
    ----------
    formula     :   str
                    Chemical formula for this species
    name        :   str
                    Common name for this species
    mu_vs_T     :   numpy array
                    Temperature array for use in viscosity interpolation
    Cp_vs_T     :   numpy array
                    Viscosity array
        
    Methods
    -------
    viscosity(self, Temperature)
        Determine the viscosity at the temperature provided to the method
    specific_heat(self, Temperature)
        Determine the viscosity at the temperature provided to the method
    plot(self, fig_num=1)
        Plots all of the data for this species on a single plot with
        several subplots
    
    '''


class H2(BaseSpecies):
    def __init__(self):
        super().__init__()
        self.name       = "Hydrogen"                        # [str] Common Name for this species
        self.formula    = "H2"                              # [str] Chemical formula for species
        self.MW         = 2.01588                           # [kg/kgmol] Molecular Mass of species
        self.MW_units   = 'kg/kgmol'                        # [str] Units for Molecuar Mass
        self.mu_units   = 'kg/m-s'                          # [str] Units for dynamic viscosity
        self.Cp_units   = 'kJ/kg-K'                         # [str] Units for Specific Heat (Cp)

        # Viscosity Data (Temperature = [K], Dynamic Viscosity = [kg/m-s])
        self.mu_vs_T = np.array([[100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 
                                  1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0],
                                 [3.779926e-6, 6.652075e-6, 1.086357e-5, 1.419705e-5, 
                                  1.708408e-5, 1.977765e-5, 2.227776e-5, 2.461417e-5, 
                                  2.684641e-5, 2.897448e-5, 3.101325e-5, 3.297763e-5, 
                                  3.489735e-5, 3.672779e-5, 3.852846e-5]])  

        # Specific Heat Data (Temperature = [K], Viscosity = [kJ/kg-K])
        self.Cp_vs_T = np.array([[300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 
                                  1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0],
                                 [14.3164, 14.4847, 14.5219, 14.5554, 14.6115, 14.7170, 
                                  14.8439, 14.6563, 15.1596, 15.3463, 15.5540, 15.7826,
                                  16.0317, 16.2456, 16.4491]])


    def viscosity(self, T):
        return np.interp(T, self.mu_vs_T[0], self.mu_vs_T[1], left=np.nan, right=np.nan)

    def specific_heat(self, T):
        return np.interp(T, self.Cp_vs_T[0], self.Cp_vs_T[1], left=np.nan, right=np.nan)


    def plot(self, fig_num=1):
        plt.figure(fig_num)
        plt.gcf().clear()
        plt.subplot(2,2,1)
        plt.plot(self.mu_vs_T[0], self.mu_vs_T[1])
        plt.plot(self.mu_vs_T[0], self.mu_vs_T[1], 'r*')
        plt.title('{0} Viscosity vs. Temperature'.format(self.name))
        plt.xlabel('Temperature [{0}]'.format('K'))
        plt.ylabel('Viscosity, $\mu$ [{}]'.format('kg/m-s'))
        plt.grid()
        
        plt.subplot(2,2,2)
        plt.plot(self.Cp_vs_T[0], self.Cp_vs_T[1])
        plt.plot(self.Cp_vs_T[0], self.Cp_vs_T[1], 'r*')
        plt.title('{0} Specific Heat vs. Temperature'.format(self.name))
        plt.xlabel('Temperature [{0}]'.format('K'))
        plt.ylabel('Specific Heat, $C_p$ [{}]'.format('kJ/kg-K'))
        plt.grid()

if __name__ == '__main__':

    fluid = H2()
    T_STP = 298.15
    
    visc = fluid.viscosity(T_STP)
    Cp = fluid.specific_heat(T_STP)
    MW = fluid.MW
    MW_units = fluid.MW_units
    mu_units = fluid.mu_units
    Cp_units = fluid.Cp_units
    
    print('Properties of the fluid {0} [{1}]'.format(fluid.name, fluid.formula))
    print('Viscosity at STP = {0:.4g} [{1}]'.format(visc, mu_units))
    print('Specific Heat at STP = {0:.4g} [{1}]'.format(Cp, Cp_units))
    print('Fluid Molecular Weight = {0:.4f} [{1}]'.format(MW, MW_units))

    fluid.plot()