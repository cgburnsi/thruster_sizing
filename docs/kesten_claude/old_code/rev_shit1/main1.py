import numpy as np
import logging
import logging.config
import json
import pathlib  # <-- Import the pathlib module
from src import utils 
from dataclasses import dataclass
from abc import ABC, abstractmethod


R_GAS = 8314.46  # Universal Gas Constant [J/(kmol*K)]



# --- Build a reliable path to the config file ---
# __file__ is the path to the current script (e.g.,.../src/main1.py)
#.parent is the directory containing the script (e.g.,.../src/)
#.parent.parent is the parent of that directory (the project root)
PROJECT_ROOT = pathlib.Path(__file__).parent.parent
CONFIG_PATH = PROJECT_ROOT / "logging_config.json"

# --- Load Logging Configuration ---
try:
    with open(CONFIG_PATH, 'rt') as f:  # <-- Use the new, full CONFIG_PATH
        config_data = json.load(f)
    logging.config.dictConfig(config_data)
except Exception as e:
    # This is the error you are seeing!
    print(f"Error loading logging config: {e}")
    # Fallback to basic config if loading fails
    logging.basicConfig(level=logging.INFO)
# ------------------------------------

# This gets the logger named "__main__"
logger = logging.getLogger(__name__)

# This gets the logger named "another.module"
other_logger = logging.getLogger("another.module")






@dataclass
class SpeciesData:
    name: str
    elemental_composition: dict
    molecular_weight: float
    
    



    
    
class AbstractThermoStrategy(ABC):
    @abstractmethod
    def get_cp(self, T: float) -> float:
        """Returns molar specific heat, Cp [J/(kmol*K)]"""
        pass
    
    @abstractmethod
    def get_h(self, T: float) -> float:
        """Returns molar enthalpy, H [J/kmol]"""
        pass
    
    @abstractmethod
    def get_s(self, T: float) -> float:
        """Returns molar entropy, S [J/(kmol*K)]"""
        pass


class NasaPolynomialStrategy(AbstractThermoStrategy):
    """
    A Concrete Strategy implementing the 7-coefficient NASA polynomial model.
    The design is modeled on implementations in libraries like pMuTT [4]
    and Cantera.[5]
    
    The polynomials are of the form [4, 6]:
    Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    S/R  = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    """
    def __init__(self, T_mid: float, a_low: np.ndarray, a_high: np.ndarray, T_low: float = 200.0, T_high: float = 6000.0):
        """
        Initializes the strategy with temperature bounds and coefficients.
        
        Args:
            T_mid (float): The midpoint temperature [K] separating low/high polys.
            a_low (np.ndarray): 7-element array for the low-T range.
            a_high (np.ndarray): 7-element array for the high-T range.
            T_low (float): Minimum valid temperature [K].
            T_high (float): Maximum valid temperature [K].
        """
        self.T_low = T_low
        self.T_mid = T_mid
        self.T_high = T_high
        self.a_low = a_low
        self.a_high = a_high

    def _get_coeffs(self, T: float) -> np.ndarray:
        """Helper method to select the correct coefficient array."""
        if T <= self.T_mid:
            return self.a_low
        else:
            return self.a_high

    def get_cp(self, T: float) -> float:
        """Returns molar specific heat, Cp [J/(kmol*K)]"""
        a = self._get_coeffs(T)
        # T_vec = [1, T, T^2, T^3, T^4]
        T_vec = np.array([T**i for i in range(5)])
        cp_r = np.dot(a[0:5], T_vec)
        return cp_r * R_GAS

    def get_h(self, T: float) -> float:
        """Returns molar enthalpy, H [J/kmol]"""
        a = self._get_coeffs(T)
        # T_vec = [1, T/2, T^2/3, T^3/4, T^4/5, 1/T]
        T_vec_terms = [T**i / (i + 1) for i in range(5)] # H/RT terms 1-5
        T_vec_terms.append(1.0 / T)                     # H/RT term 6
        T_vec = np.array(T_vec_terms)
        
        h_rt = np.dot(a[0:6], T_vec)
        return h_rt * R_GAS * T

    def get_s(self, T: float) -> float:
        """Returns molar entropy, S [J/(kmol*K)]"""
        a = self._get_coeffs(T)
        # S/R  = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
        # T_vec = [ln(T), T, T^2/2, T^3/3, T^4/4]
        T_vec_terms = [np.log(T)]
        T_vec_terms.extend([T**i / i for i in range(1, 5)])
        T_vec = np.array(T_vec_terms)
        
        s_r = np.dot(a[0:5], T_vec) + a[6]
        return s_r * R_GAS

class Species:
    """
    The Context class.
    This object holds its static data (`.data`) and its swappable
    behavior (`.thermo`). It delegates all calculation calls to the
    attached strategy object.[2]
    """
    def __init__(self, data: SpeciesData, thermo: AbstractThermoStrategy):
        self.data = data
        self.thermo = thermo

    def get_cp(self, T: float) -> float:
        """Delegates Cp calculation to its thermo strategy."""
        return self.thermo.get_cp(T)

    def get_h(self, T: float) -> float:
        """Delegates H calculation to its thermo strategy."""
        return self.thermo.get_h(T)

    def get_s(self, T: float) -> float:
        """Delegates S calculation to its thermo strategy."""
        return self.thermo.get_s(T)




class Mixture:
    """
    The Mixture (Composite/Facade) class.[7, 8]
    This class holds the system's state (T, P, Composition) and
    manages a list of 'Leaf' (Species) objects.
    It provides a high-level API for all mixture properties,
    modeling the API of professional libraries like Cantera.[9, 10]
    """
    def __init__(self, species_list):
        self.species_list = species_list
        self.n_species = len(species_list)
        
        # Create a quick lookup map for setting composition by name
        self._species_name_map = {s.data.name: i for i, s in enumerate(self.species_list)}
        
        # Initialize default state
        self._T = 298.15  # Temperature [K]
        self._P = 101325  # Pressure [Pa]
        self._X = np.zeros(self.n_species)
        if self.n_species > 0:
            self._X = 1.0  # Default to 100% of the first species

    @property
    def T(self) -> float:
        return self._T

    @T.setter
    def T(self, temperature_K: float):
        if temperature_K <= 0:
            raise ValueError("Temperature must be positive.")
        self._T = temperature_K

    @property
    def P(self) -> float:
        return self._P

    @P.setter
    def P(self, pressure_Pa: float):
        if pressure_Pa <= 0:
            raise ValueError("Pressure must be positive.")
        self._P = pressure_Pa

    @property
    def X(self) -> np.ndarray:
        """Gets the normalized mole fraction vector."""
        return self._X

    @X.setter
    def X(self, value):
        """
        Sets the mole fractions. Input can be a dict, list, or array.
        The vector is always normalized upon setting.[10]
        """
        x_vector = np.zeros(self.n_species)
        
        if isinstance(value, dict):
            for name, mole_fraction in value.items():
                if name in self._species_name_map:
                    idx = self._species_name_map[name]
                    x_vector[idx] = mole_fraction
                else:
                    # Silently ignore species not in the mixture
                    pass 
        elif isinstance(value, (list, np.ndarray)):
            if len(value)!= self.n_species:
                raise ValueError(f"Composition vector length mismatch. Expected {self.n_species}, got {len(value)}")
            x_vector = np.array(value)
        else:
            raise TypeError("Composition must be set with a dict, list, or numpy array.")

        # Normalize the vector
        total = np.sum(x_vector)
        if total == 0:
            self._X = x_vector  # Avoid division by zero
        else:
            self._X = x_vector / total

    @property
    def TPX(self):
        """Convenience getter for (T, P, X)."""
        return self._T, self._P, self._X

    @TPX.setter
    def TPX(self, values):
        """Convenience setter for (T, P, X).[10]"""
        self.T, self.P, self.X = values

    # --- Composite Property Calculations ---

    @property
    def mean_molecular_weight(self) -> float:
        """Calculates the mean molecular weight [kg/kmol]."""
        mw_vector = np.array([s.data.molecular_weight for s in self.species_list])
        return np.dot(self.X, mw_vector)

    @property
    def enthalpy_mole(self) -> float:
        """Calculates the molar enthalpy of the mixture [J/kmol]."""
        T = self.T
        h_total = 0.0
        for i, species in enumerate(self.species_list):
            # The core Composite + Strategy delegation:
            # 1. Mixture (Composite) calls Species (Leaf)
            # 2. Species (Context) delegates to Strategy
            h_species = species.get_h(T) 
            h_total += self.X[i] * h_species
        return h_total

    @property
    def cp_mole(self) -> float:
        """Calculates the molar specific heat of the mixture [J/kmol/K]."""
        T = self.T
        cp_total = 0.0
        for i, species in enumerate(self.species_list):
            cp_species = species.get_cp(T)
            cp_total += self.X[i] * cp_species
        return cp_total

    @property
    def cp_mass(self) -> float:
        """Calculates the mass specific heat of the mixture [J/kg/K]."""
        # Properties can be calculated from other properties [9, 5]
        # [J/kmol/K] / [kg/kmol] = [J/kg/K]
        return self.cp_mole / self.mean_molecular_weight


# === Part 5: The Factory Class ===

class SpeciesFactory:
    """
    The Factory pattern.[11, 12, 13]
    Decouples the client (main script) from the complex process of
    creating and assembling Species objects from their components
    (data, thermo models).
    """
    def __init__(self):
        self._species_db = {}
        # In a real application, you might call a file parser here.
        # For this example, we use a hard-coded loader.
        self._load_hardcoded_data()

    def _load_hardcoded_data(self):
        """
        This method replaces a complex file parser.
        It manually assembles and "registers" the species objects.
        """
        
        # --- N2H4 (Hydrazine) ---
        # (Dummy data for demonstration)
        n2h4_data = SpeciesData(
            name="N2H4",
            elemental_composition={"N": 2, "H": 4},
            molecular_weight=32.045  # kg/kmol
        )
        n2h4_thermo = NasaPolynomialStrategy(
            T_mid=1000.0,
            a_low=np.array([4.165, 0.015, -1.0e-5, 2.0e-9, -1.0e-13, -1000.0, 5.0]),
            a_high=np.array([5.0, 0.002, -1.0e-6, 1.0e-9, -1.0e-13, -1000.0, 5.0])
        )
        n2h4_species = Species(data=n2h4_data, thermo=n2h4_thermo)
        self._species_db["N2H4"] = n2h4_species

        # --- H2 (Hydrogen) ---
        # (Dummy data for demonstration)
        h2_data = SpeciesData(
            name="H2",
            elemental_composition={"H": 2},
            molecular_weight=2.016  # kg/kmol
        )
        h2_thermo = NasaPolynomialStrategy(
            T_mid=1000.0,
            a_low=np.array([3.5, 0.0, 0.0, 0.0, 0.0, -100.0, 2.0]),
            a_high=np.array([3.5, 0.0, 0.0, 0.0, 0.0, -100.0, 2.0])
        )
        h2_species = Species(data=h2_data, thermo=h2_thermo)
        self._species_db["H2"] = h2_species

        # --- N2 (Nitrogen) ---
        # (Dummy data for demonstration)
        n2_data = SpeciesData(
            name="N2",
            elemental_composition={"N": 2},
            molecular_weight=28.014  # kg/kmol
        )
        n2_thermo = NasaPolynomialStrategy(
            T_mid=1000.0,
            a_low=np.array([3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0]),
            a_high=np.array([3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
        )
        n2_species = Species(data=n2_data, thermo=n2_thermo)
        self._species_db["N2"] = n2_species


    def get_species(self, name: str) -> Species:
        """Retrieves a single species object from the cache."""
        species = self._species_db.get(name)
        if not species:
            raise ValueError(f"Species '{name}' not found in database.")
        return species

    def create_mixture(self, species_names) -> Mixture:
        """
        High-level factory method.
        Takes a list of species names and returns a fully-formed Mixture object.
        """
        species_objects = [self.get_species(name) for name in species_names]
        return Mixture(species_list=species_objects)







if __name__ == '__main__':
    
    # 1. Factory builds components from its data source
    print("Initializing SpeciesFactory...")
    factory = SpeciesFactory()

    # 2. Factory assembles the 'Mixture' (Composite)
    # To add a species, just add its name here (if it's in the factory)
    species_list = ['N2H4', 'H2', 'N2']
    print(f"Creating mixture with species: {species_list}")
    mixture = factory.create_mixture(species_list)

    # 3. Reactor simulation code interacts with the 'Mixture' (Facade)
    # This code is now general and species-agnostic.
    
    print("\n--- State 1 ---")
    mixture.TPX = 800.0, 5e5, {'N2H4': 1.0}

    print(f"Mixture at T={mixture.T} K, P={mixture.P:.0f} Pa")
    print(f"Composition (X): {mixture.X}")
    print(f"Mean MW (kg/kmol): {mixture.mean_molecular_weight:.4f}")
    print(f"Molar Enthalpy (J/kmol): {mixture.enthalpy_mole:,.2f}")
    print(f"Mass Cp (J/kg/K): {mixture.cp_mass:,.2f}")

    print("\n--- State 2 ---")
    # Set a new state with a different composition
    mixture.TPX = 1200.0, 5e5, {'N2H4': 0.1, 'H2': 0.5, 'N2': 0.4}
    
    print(f"Mixture at T={mixture.T} K, P={mixture.P:.0f} Pa")
    print(f"Composition (X): {mixture.X}")
    print(f"Mean MW (kg/kmol): {mixture.mean_molecular_weight:.4f}")
    print(f"Molar Enthalpy (J/kmol): {mixture.enthalpy_mole:,.2f}")
    print(f"Mass Cp (J/kg/K): {mixture.cp_mass:,.2f}")