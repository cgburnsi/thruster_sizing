# species/base_species.py

class BaseSpecies:
    def __init__(self):
        self.name = None
        self.formula = None
        self.MW = None  # [g/mol or kg/kmol]

    def viscosity(self, T):
        raise NotImplementedError("viscosity() must be implemented in the child class")

    def specific_heat(self, T):
        raise NotImplementedError("specific_heat() must be implemented in the child class")
