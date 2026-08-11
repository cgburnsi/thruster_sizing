import math

R = 8.314462618  # J/(mol·K)

class Species:
    """
    NASA-7 thermochem with multiple T ranges.
    nasa_coeffs: list of dicts like
       {"T_low": 200.0, "T_high": 1000.0, "a": [a1..a7]}
    Optionally include hf298 (J/mol) to compute absolute enthalpy.
    """
    def __init__(self, name, formula, molecular_weight, nasa_coeffs, hf298=None):
        self.name = name
        self.formula = formula
        self.molecular_weight = molecular_weight
        self.nasa_coeffs = nasa_coeffs
        self.hf298 = hf298  # standard formation enthalpy at 298.15 K, J/mol

        # Normalize/validate ordering by T_low
        self.nasa_coeffs = sorted(self.nasa_coeffs, key=lambda d: d["T_low"])

    def __repr__(self):
        return f"Species({self.name})"

    def __eq__(self, other):
        return isinstance(other, Species) and (self.name, self.formula) == (other.name, other.formula)

    def __hash__(self):
        return hash((self.name, self.formula))

    # ---- Public API ---------------------------------------------------------
    def cp(self, T):
        """Heat capacity at constant pressure [J/(mol·K)]."""
        a = self._get_coeffs_for_T(T)
        # NASA-7: Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        t, t2, t3, t4 = T, T*T, T*T*T, T*T*T*T
        cp_over_R = a[0] + a[1]*t + a[2]*t2 + a[3]*t3 + a[4]*t4
        return cp_over_R * R

    def h_sensible(self, T):
        """Sensible enthalpy relative to the polynomial reference [J/mol]."""
        a = self._get_coeffs_for_T(T)
        # NASA-7: H/RT = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
        t, t2, t3, t4 = T, T*T, T*T*T, T*T*T*T
        h_over_RT = (a[0]
                     + a[1]*t/2.0
                     + a[2]*t2/3.0
                     + a[3]*t3/4.0
                     + a[4]*t4/5.0
                     + a[5]/t)
        return h_over_RT * R * T

    def s(self, T):
        """Entropy [J/(mol·K)]."""
        a = self._get_coeffs_for_T(T)
        # NASA-7: S/R = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7
        t, t2, t3, t4 = T, T*T, T*T*T, T*T*T*T
        s_over_R = (a[0]*math.log(t)
                    + a[1]*t
                    + a[2]*t2/2.0
                    + a[3]*t3/3.0
                    + a[4]*t4/4.0
                    + a[6])
        return s_over_R * R

    def h(self, T, use_hf298=True):
        """
        Absolute enthalpy at T [J/mol].
        If hf298 is provided and use_hf298=True, returns:
            hf298 + (h_sensible(T) - h_sensible(298.15))
        Otherwise returns the polynomial's absolute h(T) (which already
        includes the a6 term matching reference conventions of that dataset).
        """
        if use_hf298 and self.hf298 is not None:
            return self.hf298 + (self.h_sensible(T) - self.h_sensible(298.15))
        # fallback: direct polynomial absolute enthalpy
        return self.h_sensible(T)

    # ---- Internals ----------------------------------------------------------
    def _get_coeffs_for_T(self, T):
        # choose interval; clamp to nearest if out of bounds
        for blk in self.nasa_coeffs:
            if blk["T_low"] <= T <= blk["T_high"]:
                return blk["a"]
        # If not found, clamp to closest block
        if T < self.nasa_coeffs[0]["T_low"]:
            return self.nasa_coeffs[0]["a"]
        return self.nasa_coeffs[-1]["a"]


class Reaction:
    """
    Stoichiometric reaction with species stoich in mol units.
    reactants/products are dicts: {Species: nu}
    """
    def __init__(self, reactants, products):
        self.reactants = dict(reactants)
        self.products = dict(products)

    def enthalpy_change(self, T, use_hf298=True):
        """
        ΔH°(T) [J/mol reaction] using absolute species enthalpies at T.
        Positive = endothermic in this sign convention.
        """
        h_prod = sum(nu * sp.h(T, use_hf298=use_hf298) for sp, nu in self.products.items())
        h_reac = sum(nu * sp.h(T, use_hf298=use_hf298) for sp, nu in self.reactants.items())
        return h_prod - h_reac


# Example (placeholder) NASA-7 blocks
H2O_blocks = [
    {"T_low": 200.0,  "T_high": 1000.0, "a": [3.38684, 3.47498e-03, -6.3547e-06, 6.96858e-09, -2.50659e-12, -3.02081e+04, 2.59023]},
    {"T_low": 1000.0, "T_high": 6000.0, "a": [2.67215, 3.05629e-03, -8.73026e-07, 1.20099e-10, -6.39162e-15, -2.98992e+04, 6.86282]}
]
# If you have ΔHf°(298.15 K) for H2O(g):
hf298_H2O = -241_826.0  # J/mol (example)

water_vapor = Species("water_vapor", "H2O", 18.01528, H2O_blocks, hf298=hf298_H2O)

# Do the same for H2 and O2 with their proper coeffs and hf298 values (elements usually 0.0 J/mol).
H2_blocks = [
    {"T_low": 200.0, "T_high": 1000.0, "a": [2.34433, 7.98052e-03, -1.94782e-05, 2.01572e-08, -7.37612e-12, -9.17935e+02, 6.83010]},
    {"T_low": 1000.0, "T_high": 6000.0, "a": [3.33728, -4.94025e-05, 4.99457e-07, -1.79566e-10, 2.00255e-14, -9.50159e+02, -3.20502]}
]
O2_blocks = [
    {"T_low": 200.0, "T_high": 1000.0, "a": [3.78246, -2.99673e-03, 9.84730e-06, -9.68130e-09, 3.24373e-12, -1.06394e+03, 3.65768]},
    {"T_low": 1000.0, "T_high": 6000.0, "a": [3.28254, 1.48309e-03, -7.57967e-07, 2.09471e-10, -2.16718e-14, -1.08846e+03, 5.45323]}
]

hydrogen = Species("hydrogen", "H2", 2.01588, H2_blocks, hf298=0.0)
oxygen   = Species("oxygen",   "O2", 31.9988, O2_blocks, hf298=0.0)

#  H2 + 0.5 O2 -> H2O(g)
h2_combustion = Reaction(reactants={hydrogen: 1.0, oxygen: 0.5},
                         products={water_vapor: 1.0})

T = 1000.0
print(f"H2O sensible h(1000 K): {water_vapor.h_sensible(T):.1f} J/mol")
print(f"H2O absolute h°(1000 K): {water_vapor.h(T):.1f} J/mol")
print(f"ΔH°_rxn(1000 K): {h2_combustion.enthalpy_change(T):.1f} J/mol")




if __name__ == '__main__':
    pass
