
class Reaction:
    def __init__(self, reactants, products, arrhenius_params):
        # reactants/products are dicts: {'N2H4': 1} -> {'NH3': 4/3, 'N2': 1/3}
        self.reactants = reactants 
        self.products = products
        self.A = arrhenius_params['A']
        self.b = arrhenius_params['b']
        self.Ea = arrhenius_params['Ea']

    def get_rate(self, temperature, concentrations):
        # Generic Arrhenius calculation: k = A * T^b * exp(-Ea/RT)
        # Rate = k * product(conc^order)
        pass
    
    def get_stoichiometry_vector(self, species_list):
        # Returns array of coefficients [-1, 0, +1.33, +0.33] matching species_list
        pass





if __name__ == '__main__':
    pass


