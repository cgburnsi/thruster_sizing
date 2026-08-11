#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
chemical_mixture.py

A dependency-free (Numpy OK) Python framework for calculating
thermodynamic and transport properties of multi-phase chemical mixtures.

Modeled after Cantera  and NASA CEA.

Usage:
    import chemical_mixture as cm
    
    mix = cm.Mixture(['N2', 'O2'])
    mix.set_state(T=500.0, P=101325.0, X={'N2': 0.79, 'O2': 0.21})
    
    print(f"Phase: {mix.phase}")
    print(f"Cp (mass): {mix.cp_mass} J/kg*K")
    print(f"Viscosity: {mix.viscosity} Pa*s")
"""

import numpy as np

# --- Section I: Constants and Solvers ---

# --- Physical Constants ---
R_GAS = 8.31446261815324  # J/(mol*K), Universal Gas Constant
P_ATM = 101325.0          # Pa, Standard Atmosphere
T_STD = 298.15            # K, Standard Temperature

class _Solvers:
    """
    Contains static methods for numerical root finding, implemented using
    only numpy to satisfy the dependency-free constraint.
    """
    @staticmethod
    def newton_1d(f, f_prime, x0, tol=1e-8, max_iter=50, bounds=None):
        """
        1D Newton-Raphson solver.
        [8, 15]
        """
        x = x0
        for _ in range(max_iter):
            fx = f(x)
            if abs(fx) < tol:
                return x
            
            fpx = f_prime(x)
            if fpx == 0:
                # Use bisection step if derivative is zero
                if bounds:
                    return _Solvers.bisection(f, bounds, bounds, tol, 20)
                return None
            
            x_new = x - fx / fpx
            
            # Simple bounding
            if bounds:
                x_new = max(bounds, min(bounds, x_new))

            if abs(x_new - x) < tol:
                return x_new
            x = x_new
            
        return None # Failed to converge

    @staticmethod
    def bisection(f, a, b, tol=1e-8, max_iter=100):
        """
        1D Bisection solver.
        
        """
        fa = f(a)
        fb = f(b)
        
        if fa * fb >= 0:
            if abs(fa) < tol: return a
            if abs(fb) < tol: return b
            return None 

        for _ in range(max_iter):
            c = (a + b) / 2.0
            fc = f(c)
            
            if abs(b - a) / 2.0 < tol or fc == 0:
                return c
            
            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
                
        return (a + b) / 2.0

    @staticmethod
    def solve_cubic_real(a, b, c, d):
        """
        Finds all real roots of a cubic polynomial ax^3 + bx^2 + cx + d = 0
        using the analytical (Cardano's) method.
        
        """
        # Normalize to depressed cubic: t^3 + p*t + q = 0
        # t = x + b/(3a)
        p = (3.0*a*c - b*b) / (3.0*a*a)
        q = (2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d) / (27.0*a*a*a)
        
        # Discriminant/4
        delta = (q / 2.0)**2 + (p / 3.0)**3
        
        roots =
        shift = b / (3.0 * a)

        if delta >= 0:
            # One real root (or 3 real, 2 equal)
            D_sqrt = np.sqrt(delta)
            A = np.cbrt(-q / 2.0 + D_sqrt)
            B = np.cbrt(-q / 2.0 - D_sqrt)
            root = A + B
            roots.append(root - shift)
            
        else:
            # Three real roots
            phi = np.arccos(-q / (2.0 * np.sqrt(-(p/3.0)**3)))
            r_term = 2.0 * np.sqrt(-p / 3.0)
            r1 = r_term * np.cos(phi / 3.0)
            r2 = r_term * np.cos((phi + 2.0 * np.pi) / 3.0)
            r3 = r_term * np.cos((phi + 4.0 * np.pi) / 3.0)
            roots = [r1 - shift, r2 - shift, r3 - shift]

        return sorted([r for r in roots if np.isreal(r)])


# --- Section II: Species Database and Class ---

# Atomic diffusion volumes (cm^3/mol) for Fuller-Schettler-Giddings
# [23, 24]
_FSG_ATOMIC_VOLUMES = {
    'C': 15.9, 'H': 2.31, 'O': 6.11, 'N': 4.54, 'Ar': 16.2, 'He': 2.67,
    # Simple molecules
    'H2': 6.12, 'N2': 18.5, 'O2': 16.3, 'CO2': 26.9, 'H2O': 13.1, 'CH4': 24.42
}

# Master database of all species parameters.
# This serves the role of Cantera's.yaml or CHEMKIN's.dat files [16]
SPECIES_DATABASE = {
    'N2': {
        'MW': 28.0134, # g/mol
        'Tc': 126.2,   # K
        'Pc': 3.39e6,  # Pa
        'omega': 0.038,
        'Tb': 77.36,   # K
        'fsg_volumes': {'N2': 1},
        'thermo_coeffs':
            [3.531005280E+00, -1.236609870E-04, -5.029994330E-08, 1.934058680E-10, -1.603875380E-13, -1.041657350E+03, -3.213437150E+00],
            # High-T: 1000-5000 K
            [2.952576370E+00, 1.396669450E-03, -4.926285140E-07, 1.088420030E-10, -9.600700400E-15, -9.227876100E+02, 3.987742960E+00],
        'gas_visc_coeffs': [-11.45, 0.9023, -0.0159], # ln(mu) = A + B*ln(T) + C*(ln(T))^2
        'gas_k_coeffs': [-4.897, 1.139, -0.024],
        'liq_visc_coeffs': [1.97e-5, 332.0] # ln(mu_L) = ln(A) + B/T 
    },
    'O2': {
        'MW': 31.9988,
        'Tc': 154.6,
        'Pc': 5.043e6,
        'omega': 0.022,
        'Tb': 90.18,
        'fsg_volumes': {'O2': 1},
        'thermo_coeffs':,
            [3.298677000E+00, 1.315136350E-03, -1.868777900E-07, -9.817088470E-11, 3.013726560E-15, -1.018318040E+03, 3.486252990E+00],
        'gas_visc_coeffs': [-11.16, 0.9096, -0.0177],
        'gas_k_coeffs': [-5.485, 1.325, -0.043],
        'liq_visc_coeffs': [1.60e-5, 309.0]
    },
    'H2O': {
        'MW': 18.01528,
        'Tc': 647.1,
        'Pc': 22.064e6,
        'omega': 0.344,
        'Tb': 373.15,
        'fsg_volumes': {'H2O': 1},
        'thermo_coeffs':,
            [4.198640560E+00, 2.306806040E-04, 5.768239030E-07, -5.510032510E-10, 1.802106950E-13, -3.020811320E+04, -8.490322080E-01],
        'gas_visc_coeffs': [-13.43, 1.371, -0.038],
        'gas_k_coeffs': [-6.621, 1.764, -0.053],
        'liq_visc_coeffs': [-11.69, 2106.0] # ln(mu_L) = A + B/T
    },
    'CO2': {
        'MW': 44.01,
        'Tc': 304.13,
        'Pc': 7.377e6,
        'omega': 0.224,
        'Tb': 194.7, # Sublimation
        'fsg_volumes': {'CO2': 1},
        'thermo_coeffs':,
            [4.632313180E+00, 2.744158040E-03, -1.218336340E-06, 2.822701830E-10, -2.627672100E-14, -4.891907680E+04, -1.531718010E+00],
        'gas_visc_coeffs': [-13.51, 1.341, -0.040],
        'gas_k_coeffs': [-8.421, 1.954, -0.076],
        'liq_visc_coeffs': [-8.70, 775.0]
    },
    'Ar': {
        'MW': 39.948,
        'Tc': 150.86,
        'Pc': 4.898e6,
        'omega': -0.001,
        'Tb': 87.3,
        'fsg_volumes': {'Ar': 1},
        'thermo_coeffs':,
            [2.500000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, 0.000000000E+00, -7.453750000E+02, 4.379674910E+00],
        'gas_visc_coeffs': [-11.02, 0.8615, -0.015],
        'gas_k_coeffs': [-5.454, 0.811, -0.013],
        'liq_visc_coeffs': [-9.60, 237.0]
    }
}


class Species:
    """
    Data container for a single species, populated from the SPECIES_DATABASE.
    """
    def __init__(self, name):
        if name not in SPECIES_DATABASE:
            raise ValueError("Species '{}' not found in database.".format(name))
            
        data = SPECIES_DATABASE[name]
        
        self.name = name
        self.mw = data
        self.tc = data
        self.pc = data['Pc']
        self.omega = data['omega']
        self.tb = data
        self.fsg_volumes = data['fsg_volumes']
        
        # NASA 7-term polynomial coefficients [25]
        self.thermo_t_mid = 1000.0 # Standard midpoint
        self.thermo_coeffs_low = data['thermo_coeffs']
        self.thermo_coeffs_high = data['thermo_coeffs']
        
        self.gas_visc_coeffs = data['gas_visc_coeffs']
        self.gas_k_coeffs = data['gas_k_coeffs']
        self.liq_visc_coeffs = data['liq_visc_coeffs']
        
        # Pre-calculate FSG diffusion volume
        self.fsg_v = self._calculate_fsg_volume()

    def _calculate_fsg_volume(self):
        """
        Calculates the total diffusion volume for the molecule.
        [23, 26]
        """
        vol = 0.0
        for atom, count in self.fsg_volumes.items():
            if atom in _FSG_ATOMIC_VOLUMES:
                vol += _FSG_ATOMIC_VOLUMES[atom] * count
            else:
                # This handles cases where the key is the molecule itself, e.g. 'N2'
                if self.name in _FSG_ATOMIC_VOLUMES:
                    return _FSG_ATOMIC_VOLUMES[self.name]
                else:
                    raise Warning("Warning: Unknown atom '{}' in FSG volumes for {}".format(atom, self.name))
        return vol


# --- Helper Data Classes ---

class MixtureThermoProperties:
    """Data container for mixture thermodynamic properties."""
    def __init__(self):
        self.mw_mix = 0.0   # g/mol
        self.rho_mass = 0.0 # kg/m^3
        self.rho_mole = 0.0 # mol/m^3
        self.cp_mole = 0.0  # J/(mol*K)
        self.cp_mass = 0.0  # J/(kg*K)
        self.h_mole = 0.0   # J/mol
        self.h_mass = 0.0   # J/kg
        self.s_mole = 0.0   # J/(mol*K)
        self.s_mass = 0.0   # J/(kg*K)

class MixtureTransportProperties:
    """Data container for mixture transport properties."""
    def __init__(self):
        self.viscosity = 0.0            # Pa*s
        self.thermal_conductivity = 0.0 # W/(m*K)
        self.mix_diff_coeffs = None     # m^2/s (vector)

class PhaseProperties:
    """Data container for properties of a single phase (L or V)."""
    def __init__(self):
        self.composition = None
        self.thermo = MixtureThermoProperties()
        self.transport = MixtureTransportProperties()


# --- Section III & VI: Ideal Gas Phase Model (Thermo & Transport) ---

class IdealGasPhase:
    """
    Calculates thermodynamic and transport properties for an ideal gas
    mixture.
    """
    def __init__(self):
        pass # No state stored

    # --- Thermo Methods (Section III) ---

    def _get_species_thermo(self, species, T):
        """
        Calculates dimensionless thermo properties for one species
        using the NASA 7-term polynomial format. 
        """
        if T > species.thermo_t_mid:
            coeffs = species.thermo_coeffs_high
        else:
            coeffs = species.thermo_coeffs_low
            
        a1, a2, a3, a4, a5, a6, a7 = coeffs
        
        T_vec = np.array()
        
        # Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        cp_R = np.dot(T_vec, np.array([a1, a2, a3, a4, a5]))
        
        # H/RT = a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
        h_RT = (a1 + a2*T/2.0 + a3*(T**2)/3.0 + a4*(T**3)/4.0 + 
                a5*(T**4)/5.0 + a6/T)
        
        # S/R = a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
        s_R = (a1*np.log(T) + a2*T + a3*(T**2)/2.0 + a4*(T**3)/3.0 + 
               a5*(T**4)/4.0 + a7)

        return cp_R, h_RT, s_R

    def calculate_thermo_properties(self, T, P, X, species_data):
        """
        Calculates the mixture thermodynamic properties for the ideal gas.
        Returns a MixtureThermoProperties data object.
        """
        n_species = len(species_data)
        
        cp_R_vec = np.zeros(n_species)
        h_RT_vec = np.zeros(n_species)
        s_R_vec = np.zeros(n_species)
        mw_vec = np.zeros(n_species)

        for i in range(n_species):
            species = species_data[i]
            mw_vec[i] = species.mw
            cp_R, h_RT, s_R = self._get_species_thermo(species, T)
            cp_R_vec[i] = cp_R
            h_RT_vec[i] = h_RT
            s_R_vec[i] = s_R
        
        # Mixture average molecular weight
        mw_mix = np.dot(X, mw_vec) # g/mol
        mw_mix_kg = mw_mix / 1000.0 # kg/mol
        
        # Molar properties
        cp_mole = R_GAS * np.dot(X, cp_R_vec)  # J/(mol*K)
        h_mole = R_GAS * T * np.dot(X, h_RT_vec) # J/mol
        
        # Add entropy of mixing
        s_mole_ideal = R_GAS * np.dot(X, s_R_vec)
        s_mix = -R_GAS * np.dot(X[X > 0], np.log(X[X > 0])) # Avoid log(0)
        s_mole = s_mole_ideal + s_mix # J/(mol*K)
        
        # Mass-specific properties
        cp_mass = cp_mole / mw_mix_kg # J/(kg*K)
        h_mass = h_mole / mw_mix_kg   # J/kg
        s_mass = s_mole / mw_mix_kg   # J/(kg*K)
        
        # Density (ideal gas law)
        rho_mole = P / (R_GAS * T) # mol/m^3
        rho_mass = rho_mole * mw_mix_kg # kg/m^3
        
        props = MixtureThermoProperties()
        props.mw_mix = mw_mix
        props.rho_mass = rho_mass
        props.rho_mole = rho_mole
        props.cp_mole = cp_mole
        props.h_mole = h_mole
        props.s_mole = s_mole
        props.cp_mass = cp_mass
        props.h_mass = h_mass
        props.s_mass = s_mass
        
        return props

    # --- Transport Methods (Section VI) ---

    def _get_species_gas_transport(self, species, T):
        """
        Calculates pure-species gas viscosity and thermal conductivity
        from polynomial fits. 
        """
        logT = np.log(T)
        
        # ln(mu) = c0 + c1*ln(T) + c2*(ln(T))^2
        c_visc = species.gas_visc_coeffs
        log_visc = c_visc + c_visc * logT + c_visc * logT**2
        mu = np.exp(log_visc) # Pa*s
        
        # ln(k) = c0 + c1*ln(T) + c2*(ln(T))^2
        c_k = species.gas_k_coeffs
        log_k = c_k + c_k * logT + c_k * logT**2
        k = np.exp(log_k) # W/(m*K)
        
        return mu, k

    def _calculate_wilke_phi_matrix(self, mus, mws):
        """
        Calculates the NxN phi_ij matrix for Wilke/Mason-Saxena.
        [42, 43]
        """
        n = len(mus)
        # Use numpy broadcasting to build matrices
        # mu_i / mu_j
        mu_i_over_j = mus[:, np.newaxis] / (mus[np.newaxis, :] + 1e-100)
        # M_j / M_i
        mw_j_over_i = mws[np.newaxis, :] / (mws[:, np.newaxis] + 1e-100)
        # M_i / M_j
        mw_i_over_j = mws[:, np.newaxis] / (mws[np.newaxis, :] + 1e-100)
        
        # Numerator: [1 + (mu_i/mu_j)^0.5 * (M_j/M_i)^0.25]^2
        numerator = (1.0 + np.sqrt(mu_i_over_j) * (mw_j_over_i)**0.25)**2
        
        # Denominator: [8 * (1 + M_i/M_j)]^0.5
        denominator = np.sqrt(8.0 * (1.0 + mw_i_over_j))
        
        phi = numerator / (denominator + 1e-100)
        return phi

    def _get_binary_diffs(self, T, P, species_data):
        """
        Calculates the NxN matrix of binary diffusion coefficients (D_ij)
        using the Fuller-Schettler-Giddings (FSG) equation.
        
        """
        n = len(species_data)
        
        mws = np.array([s.mw for s in species_data])
        vols = np.array([s.fsg_v for s in species_data])
        
        # Pre-calculate T^1.75 and P (in Pa)
        T_pow = T**1.75
        P_pa = P
        
        # Use broadcasting
        # sqrt(1/M_i + 1/M_j)
        mw_term = np.sqrt(1.0/mws[:, np.newaxis] + 1.0/mws[np.newaxis, :])
        # [(V_i)^1/3 + (V_j)^1/3]^2
        vol_term = (vols[:, np.newaxis]**(1.0/3.0) + vols[np.newaxis, :]**(1.0/3.0))**2
        
        # D_ij in m^2/s
        D_ij = (1.0e-7 * T_pow * mw_term) / (P_pa * vol_term + 1e-100)
        
        return D_ij

    def calculate_transport_properties(self, T, P, X, species_data):
        """
        Calculates the mixture-averaged transport properties for the ideal gas.
        Returns a MixtureTransportProperties data object.
        """
        n_species = len(species_data)
        
        mu_vec = np.zeros(n_species)
        k_vec = np.zeros(n_species)
        mw_vec = np.array([s.mw for s in species_data])

        for i in range(n_species):
            mu_i, k_i = self._get_species_gas_transport(species_data[i], T)
            mu_vec[i] = mu_i
            k_vec[i] = k_i
            
        # --- Viscosity (Wilke) & Thermal Conductivity (Mason-Saxena) ---
        phi_matrix = self._calculate_wilke_phi_matrix(mu_vec, mw_vec)
        
        # Denominator for Wilke/Mason-Saxena: SUM(j) [X_j * phi_ij]
        phi_sum_vec = np.dot(phi_matrix, X)
        
        # mu_m = SUM(i) [X_i * mu_i / (phi_sum_vec)_i]
        mu_mix = np.sum(X * mu_vec / (phi_sum_vec + 1e-100))
        
        # k_m = SUM(i) [X_i * k_i / (phi_sum_vec)_i]
        k_mix = np.sum(X * k_vec / (phi_sum_vec + 1e-100))
        
        # --- Diffusion Coefficients (FSG + Mixture-Averaged) ---
        D_ij_matrix = self._get_binary_diffs(T, P, species_data)
        
        D_mix_vec = np.zeros(n_species)
        for k in range(n_species):
            # Denominator: SUM(j!=k)
            sum_term = 0.0
            for j in range(n_species):
                if j == k:
                    continue
                sum_term += X[j] / (D_ij_matrix[k, j] + 1e-100)
                
            if sum_term < 1e-100:
                # Handle pure species or binary case edge
                D_mix_vec[k] = 0.0
            else:
                # D_k_m = (1 - X_k) / SUM(j!=k)
                D_mix_vec[k] = (1.0 - X[k]) / sum_term
                
        props = MixtureTransportProperties()
        props.viscosity = mu_mix
        props.thermal_conductivity = k_mix
        props.mix_diff_coeffs = D_mix_vec
        
        return props


# --- Section IV & VII: Real Fluid (PR) Phase (VLE & Liquid Transport) ---

class PengRobinsonPhase:
    """
    Implements the Peng-Robinson (PR) Equation of State.
    This class is used for two main purposes:
    1. VLE: Calculating fugacity coefficients for liquid and vapor phases.
    2. Liquid Transport: Calculating liquid-phase transport properties.
    """
    def __init__(self):
        self.R = R_GAS

    # --- VLE Methods (Section IV) ---

    def _calculate_mixture_params(self, T, X, species_data):
        """
        Calculates PR EoS mixture parameters a_m and b_m.
        [18, 32]
        """
        n_species = len(species_data)
        a_i_vec = np.zeros(n_species)
        b_i_vec = np.zeros(n_species)
        
        for i in range(n_species):
            s = species_data[i]
            Tr = T / s.tc
            
            # Pure component parameters
            a_c = 0.45724 * (self.R * s.tc)**2 / s.pc
            b_c = 0.07780 * (self.R * s.tc) / s.pc
            kappa = 0.37464 + 1.54226 * s.omega - 0.26992 * s.omega**2
            alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
            
            a_i_vec[i] = a_c * alpha
            b_i_vec[i] = b_c
        
        # van der Waals mixing rules [32]
        # Assume k_ij = 0 
        # a_m = SUM(i, j) [X_i * X_j * sqrt(a_i * a_j)]
        a_ij_sqrt = np.sqrt(a_i_vec)
        a_m = (np.dot(X, a_ij_sqrt))**2 
        
        # b_m = SUM(i) [X_i * b_i]
        b_m = np.dot(X, b_i_vec)
        
        return a_m, b_m, a_i_vec, b_i_vec

    def _solve_Z(self, T, P, a_m, b_m):
        """
        Solves the PR EoS cubic polynomial for Z-factor.
        Returns a sorted list of real roots.
        """
        A = a_m * P / (self.R * T)**2
        B = b_m * P / (self.R * T)
        
        # Z^3 - (1-B)Z^2 + (A-2B-3B^2)Z - (AB-B^2-B^3) = 0
        c2 = -(1.0 - B)
        c1 = (A - 2.0*B - 3.0*B**2)
        c0 = -(A*B - B**2 - B**3)
        
        roots = _Solvers.solve_cubic_real(1.0, c2, c1, c0)
        return roots

    def calculate_fugacity_coeffs(self, T, P, X, species_data):
        """
        Calculates the fugacity coefficients for both liquid (min Z)
        and vapor (max Z) phases. 
        """
        a_m, b_m, a_i_vec, b_i_vec = self._calculate_mixture_params(T, X, species_data)
        Z_roots = self._solve_Z(T, P, a_m, b_m)
        
        if not Z_roots:
            return None, None

        Z_L = Z_roots
        Z_V = Z_roots[-1]
        
        phi_L = self._calc_phi_vector(T, P, X, species_data, Z_L, a_m, b_m, a_i_vec, b_i_vec)
        phi_V = self._calc_phi_vector(T, P, X, species_data, Z_V, a_m, b_m, a_i_vec, b_i_vec)
        
        return phi_L, phi_V

    def _calc_phi_vector(self, T, P, X, species_data, Z, a_m, b_m, a_i_vec, b_i_vec):
        """
        Helper function to compute the vector of fugacity coefficients
        for a given phase (Z-factor).
        """
        n_species = len(species_data)
        ln_phi_vec = np.zeros(n_species)
        
        A = a_m * P / (self.R * T)**2
        B = b_m * P / (self.R * T)

        # Pre-calculate cross-term: (2 * SUM(j) [X_j * a_ij]) / a_m
        # With k_ij=0, a_ij = sqrt(a_i * a_j)
        a_i_sqrt = np.sqrt(a_i_vec)
        sum_term = np.dot(X, a_i_sqrt)
        cross_term_vec = (2.0 * a_i_sqrt * sum_term) / (a_m + 1e-100)

        log_term_arg = (Z + (1.0 + np.sqrt(2.0)) * B) / (Z + (1.0 - np.sqrt(2.0)) * B)
        # Avoid log(negative) if Z/B are pathological
        if log_term_arg < 1e-100: log_term_arg = 1e-100
        
        log_term = np.log(log_term_arg)
        
        C1 = (Z - 1.0)
        C2_arg = Z - B
        if C2_arg < 1e-100: C2_arg = 1e-100
        C2 = np.log(C2_arg)
        C3 = A / (2.0 * np.sqrt(2.0) * B + 1e-100)
        
        for i in range(n_species):
            bi_b_ratio = (b_i_vec[i] / (b_m + 1e-100))
            term1 = bi_b_ratio * C1
            term2 = C2
            term3 = C3 * (cross_term_vec[i] - bi_b_ratio) * log_term
            
            ln_phi_vec[i] = term1 - term2 - term3
            
        return np.exp(ln_phi_vec)

    # --- Liquid Transport Methods (Section VII) ---

    def _get_pure_liquid_transport(self, T, species):
        """
        Calculates pure-species liquid viscosity and thermal conductivity.
        """
        # Viscosity: Arrhenius-type ln(mu_L) = A + B/T 
        c_visc = species.liq_visc_coeffs
        mu_L = np.exp(c_visc + c_visc / T)
        
        # Thermal Conductivity: Sato-Riedel 
        Tr = T / species.tc
        Tbr = species.tb / species.tc
        if Tr >= 1.0: # Above critical temp
            k_L = 1e-10 # Return small number
        else:
            term_num = 3.0 + 20.0 * (1.0 - Tr)**(2.0/3.0)
            term_den = 3.0 + 20.0 * (1.0 - Tbr)**(2.0/3.0)
            k_L = (1.1053 / np.sqrt(species.mw)) * (term_num / term_den)
            if k_L < 0: k_L = 1e-10
            
        return mu_L, k_L

    def calculate_liquid_transport(self, T, P, X, species_data):
        """
        Calculates liquid mixture transport properties.
        """
        n_species = len(species_data)
        mu_L_vec = np.zeros(n_species)
        k_L_vec = np.zeros(n_species)
        mw_vec = np.array([s.mw for s in species_data])
        
        for i in range(n_species):
            mu_L_i, k_L_i = self._get_pure_liquid_transport(T, species_data[i])
            mu_L_vec[i] = mu_L_i
            k_L_vec[i] = k_L_i
            
        # --- Viscosity: Logarithmic mixing by mass fraction  ---
        mw_mix = np.dot(X, mw_vec)
        if mw_mix == 0:
            W = X # Failsafe
        else:
            W = (X * mw_vec) / mw_mix # Mass fractions
        
        # ln(mu_mix) = SUM(w_i * ln(mu_i))
        mu_mix = np.exp(np.dot(W, np.log(mu_L_vec + 1e-100)))
        
        # --- Thermal Conductivity: Ideal (linear) mixing by mole fraction [52] ---
        # k_mix = SUM(x_i * k_i)
        k_mix = np.dot(X, k_L_vec)
        
        props = MixtureTransportProperties()
        props.viscosity = mu_mix
        props.thermal_conductivity = k_mix
        props.mix_diff_coeffs = np.zeros(n_species) # Diffusion in liquids not implemented
        
        return props


# --- Section V & VIII: Main Mixture Class (VLE & Facade) ---

class Mixture:
    """
    Main user-facing class for calculating mixture properties.
    Modeled after the Cantera Solution object.
    """
    def __init__(self, species_names):
        """
        Initializes the mixture with a list of species names.
        """
        self.species_names = species_names
        self.species_data =
        self.n_species = len(species_names)
        
        # Instantiate the phase models
        self.gas_phase_model = IdealGasPhase()
        self.pr_phase_model = PengRobinsonPhase()
        
        # State variables
        self.T = None
        self.P = None
        self.X = None # Overall (feed) composition
        
        # Phase state results
        self.phase = 'unknown' # 'gas', 'liquid', 'VLE'
        self.vapor_fraction = 0.0 # Phi_v
        
        # Phase-specific properties
        self.liquid_properties = PhaseProperties()
        self.vapor_properties = PhaseProperties()
        
        # Overall mixture properties (averaged)
        self.thermo = MixtureThermoProperties()
        self.transport = MixtureTransportProperties()

    def _x_dict_to_vec(self, X_dict):
        """Converts a composition dict to a numpy vector."""
        X_vec = np.zeros(self.n_species)
        for i, name in enumerate(self.species_names):
            X_vec[i] = X_dict.get(name, 0.0)
        
        # Normalize
        s = np.sum(X_vec)
        if s == 0:
            raise ValueError("Composition cannot be all zero.")
        X_vec /= s
        return X_vec

    # --- Section V: VLE Flash Logic ---
    
    def _rachford_rice_obj(self, phi_v, Z, K):
        """
        Rachford-Rice objective function f(phi_v).
        
        """
        return np.sum(Z * (K - 1.0) / (1.0 + phi_v * (K - 1.0)))
    
    def _rachford_rice_deriv(self, phi_v, Z, K):
        """
        Derivative of the Rachford-Rice objective function f'(phi_v).
        
        """
        return np.sum(-Z * (K - 1.0)**2 / (1.0 + phi_v * (K - 1.0))**2)

    def _wilson_k_guess(self, T, P):
        """
        Generates an initial guess for K-values using Wilson's correlation.
        """
        K = np.zeros(self.n_species)
        for i, s in enumerate(self.species_data):
            Tr = T / s.tc
            K[i] = (s.pc / P) * np.exp(5.37 * (1.0 + s.omega) * (1.0 - 1.0 / Tr))
        return K
        
    def _perform_tp_flash(self, T, P, Z):
        """
        Performs the main T-P flash calculation to find the phase state.
        
        """
        # 1. Initial K-value guess
        K = self._wilson_k_guess(T, P)
        
        # Check phase bounds
        f0 = self._rachford_rice_obj(0.0, Z, K)
        f1 = self._rachford_rice_obj(1.0, Z, K)
        
        if f0 <= 0:
            return 'liquid', 0.0, Z, None # Subcooled liquid
        if f1 >= 0:
            return 'gas', 1.0, None, Z # Superheated vapor
            
        phi_v_guess = 0.5
        
        # 2. Outer Loop: Successive Substitution
        for _ in range(20): # Max VLE iterations
            
            # 3. Inner Loop: Solve Rachford-Rice for phi_v
            f_obj = lambda phi: self._rachford_rice_obj(phi, Z, K)
            f_deriv = lambda phi: self._rachford_rice_deriv(phi, Z, K)
            
            # Use Newton, fallback to bisection
            phi_v = _Solvers.newton_1d(f_obj, f_deriv, phi_v_guess, tol=1e-8, bounds=(0.0, 1.0))
            if phi_v is None:
                phi_v = _Solvers.bisection(f_obj, 0.0, 1.0, tol=1e-8)
            
            if phi_v is None:
                # Failed to solve flash, unknown state
                return 'unknown', 0.0, None, None
            
            # 4. Update phase compositions
            x_liq = Z / (1.0 + phi_v * (K - 1.0))
            x_liq /= np.sum(x_liq) # Renormalize
            
            y_vap = K * x_liq
            y_vap /= np.sum(y_vap) # Renormalize
            
            # 5. Update Fugacities (EoS Coupling)
            phi_L, phi_V = self.pr_phase_model.calculate_fugacity_coeffs(T, P, x_liq, self.species_data)
            _, phi_V_y = self.pr_phase_model.calculate_fugacity_coeffs(T, P, y_vap, self.species_data)
            
            if phi_L is None or phi_V_y is None:
                 return 'unknown', 0.0, None, None
            
            # 6. Update K-values
            K_new = phi_L / phi_V_y
            
            # 7. Check Convergence
            err = np.sum((np.log(K_new + 1e-100) - np.log(K + 1e-100))**2)
            if err < 1e-10:
                return 'VLE', phi_v, x_liq, y_vap
                
            K = K_new
            phi_v_guess = phi_v

        return 'unknown', 0.0, None, None # Failed to converge VLE

    # --- Main Public API ---

    def set_state(self, T, P, X):
        """
        Sets the thermodynamic state of the mixture.
        This is the main calculation method.
        """
        self.T = T
        self.P = P
        self.X = self._x_dict_to_vec(X)
        
        # 1. Determine phase state
        phase, phi_v, x_liq, y_vap = self._perform_tp_flash(T, P, self.X)
        self.phase = phase
        self.vapor_fraction = phi_v
        
        # 2. Calculate properties based on phase
        if phase == 'gas':
            self.vapor_fraction = 1.0
            self.vapor_properties.composition = self.X
            self.vapor_properties.thermo = self.gas_phase_model.calculate_thermo_properties(T, P, self.X, self.species_data)
            self.vapor_properties.transport = self.gas_phase_model.calculate_transport_properties(T, P, self.X, self.species_data)
            
            # Set top-level properties
            self.thermo = self.vapor_properties.thermo
            self.transport = self.vapor_properties.transport
            
        elif phase == 'liquid':
            self.vapor_fraction = 0.0
            self.liquid_properties.composition = self.X
            # Use Ideal Gas thermo (standard simplification)
            self.liquid_properties.thermo = self.gas_phase_model.calculate_thermo_properties(T, P, self.X, self.species_data)
            # Use liquid transport models
            self.liquid_properties.transport = self.pr_phase_model.calculate_liquid_transport(T, P, self.X, self.species_data)

            # Set top-level properties
            self.thermo = self.liquid_properties.thermo
            self.transport = self.liquid_properties.transport
        
        elif phase == 'VLE':
            # Calculate properties for each phase
            self.vapor_properties.composition = y_vap
            self.vapor_properties.thermo = self.gas_phase_model.calculate_thermo_properties(T, P, y_vap, self.species_data)
            self.vapor_properties.transport = self.gas_phase_model.calculate_transport_properties(T, P, y_vap, self.species_data)

            self.liquid_properties.composition = x_liq
            self.liquid_properties.thermo = self.gas_phase_model.calculate_thermo_properties(T, P, x_liq, self.species_data)
            self.liquid_properties.transport = self.pr_phase_model.calculate_liquid_transport(T, P, x_liq, self.species_data)
            
            # Calculate mixture-averaged properties
            self.thermo = self._average_vle_thermo(phi_v, self.liquid_properties.thermo, self.vapor_properties.thermo)
            # Mixture-averaged transport is ill-defined and not calculated
            self.transport = MixtureTransportProperties() 
            
        else:
            raise RuntimeError("Phase calculation failed to converge.")

    def _average_vle_thermo(self, phi_v, liq_thermo, vap_thermo):
        """
        Calculates mole-fraction-weighted average properties for a VLE mixture.
        """
        phi_l = 1.0 - phi_v
        
        # Get mole-basis properties
        h_liq = liq_thermo.h_mole
        h_vap = vap_thermo.h_mole
        s_liq = liq_thermo.s_mole
        s_vap = vap_thermo.s_mole
        cp_liq = liq_thermo.cp_mole
        cp_vap = vap_thermo.cp_mole
        
        # Average on mole basis
        h_mix_mole = phi_l * h_liq + phi_v * h_vap
        s_mix_mole = phi_l * s_liq + phi_v * s_vap
        cp_mix_mole = phi_l * cp_liq + phi_v * cp_vap
        
        # Mixture MW
        mw_mix = phi_l * liq_thermo.mw_mix + phi_v * vap_thermo.mw_mix
        mw_mix_kg = mw_mix / 1000.0
        
        # Convert to mass basis
        props = MixtureThermoProperties()
        props.mw_mix = mw_mix
        props.h_mole = h_mix_mole
        props.s_mole = s_mix_mole
        props.cp_mole = cp_mix_mole
        props.h_mass = h_mix_mole / mw_mix_kg
        props.s_mass = s_mix_mole / mw_mix_kg
        props.cp_mass = cp_mix_mole / mw_mix_kg
        
        # Average density
        rho_L = liq_thermo.rho_mass
        rho_V = vap_thermo.rho_mass
        props.rho_mass = (1.0 / ( (1-phi_v)/rho_L + phi_v/rho_V ))
        
        return props


    # --- Convenience Properties (Facade)  ---
    @property
    def cp_mass(self):
        return self.thermo.cp_mass

    @property
    def cp_mole(self):
        return self.thermo.cp_mole

    @property
    def h_mass(self):
        return self.thermo.h_mass
        
    @property
    def h_mole(self):
        return self.thermo.h_mole

    @property
    def s_mass(self):
        return self.thermo.s_mass
        
    @property
    def s_mole(self):
        return self.thermo.s_mole

    @property
    def density(self):
        return self.thermo.rho_mass
        
    @property
    def mean_molecular_weight(self):
        return self.thermo.mw_mix

    @property
    def viscosity(self):
        return self.transport.viscosity

    @property
    def thermal_conductivity(self):
        return self.transport.thermal_conductivity

    @property
    def mix_diff_coeffs(self):
        return self.transport.mix_diff_coeffs


# --- Section VIII: Usage Examples ---

if __name__ == "__main__":
    
    print("--- Example 1: Single-Phase Gas (Air) ---")
    mix_air = Mixture(['N2', 'O2', 'Ar'])
    X_air = {'N2': 0.78, 'O2': 0.21, 'Ar': 0.01}
    T_gas = 800.0 # K
    P_gas = 2.0 * P_ATM # Pa
    
    mix_air.set_state(T=T_gas, P=P_gas, X=X_air)
    
    print("State (T={:.1f} K, P={:.2e} Pa):".format(T_gas, P_gas))
    print("  Phase: {}".format(mix_air.phase))
    print("  Density: {:.4f} kg/m^3".format(mix_air.density))
    print("  Cp (mass): {:.2f} J/kg*K".format(mix_air.cp_mass))
    print("  Viscosity: {:.2e} Pa*s".format(mix_air.viscosity))
    print("  Conductivity: {:.2e} W/m*K".format(mix_air.thermal_conductivity))
    print("  D_mix (N2): {:.2e} m^2/s".format(mix_air.mix_diff_coeffs))
    print("  D_mix (O2): {:.2e} m^2/s".format(mix_air.mix_diff_coeffs))


    print("\n--- Example 2: Single-Phase Liquid (Water) ---")
    mix_water = Mixture(['H2O'])
    X_water = {'H2O': 1.0}
    T_liq = 300.0 # K
    P_liq = 1.0 * P_ATM # Pa
    
    mix_water.set_state(T=T_liq, P=P_liq, X=X_water)

    print("State (T={:.1f} K, P={:.2e} Pa):".format(T_liq, P_liq))
    print("  Phase: {}".format(mix_water.phase))
    # Note: Ideal gas density is used for thermo, not liquid density
    print("  Density (ideal gas ref): {:.4f} kg/m^3".format(mix_water.density))
    print("  Cp (mass): {:.2f} J/kg*K".format(mix_water.cp_mass))
    print("  Viscosity (liquid): {:.2e} Pa*s".format(mix_water.viscosity))
    print("  Conductivity (liquid): {:.2e} W/m*K".format(mix_water.thermal_conductivity))


    print("\n--- Example 3: Two-Phase VLE (Water + Nitrogen) ---")
    mix_vle = Mixture(['H2O', 'N2'])
    X_vle = {'H2O': 0.5, 'N2': 0.5}
    T_vle = 400.0 # K
    P_vle = 50.0e5 # 50 bar
    
    mix_vle.set_state(T=T_vle, P=P_vle, X=X_vle)
    
    print("State (T={:.1f} K, P={:.2e} Pa):".format(T_vle, P_vle))
    print("  Phase: {}".format(mix_vle.phase))
    print("  Vapor Fraction (mol): {:.4f}".format(mix_vle.vapor_fraction))
    
    # Liquid Phase
    liq_X = mix_vle.liquid_properties.composition
    print("  Liquid Phase ({:.2f} %):".format((1.0 - mix_vle.vapor_fraction)*100))
    print("    X_H2O: {:.4f}, X_N2: {:.4f}".format(liq_X, liq_X))
    print("    Viscosity: {:.2e} Pa*s".format(mix_vle.liquid_properties.transport.viscosity))

    # Vapor Phase
    vap_X = mix_vle.vapor_properties.composition
    print("  Vapor Phase ({:.2f} %):".format(mix_vle.vapor_fraction*100))
    print("    Y_H2O: {:.4f}, Y_N2: {:.4f}".format(vap_X, vap_X))
    print("    Viscosity: {:.2e} Pa*s".format(mix_vle.vapor_properties.transport.viscosity))
    
    # Overall Mixture
    print("  Overall Mixture Properties:")
    print("    h (mass): {:.2e} J/kg".format(mix_vle.h_mass))
    print("    Cp (mass): {:.2f} J/kg*K".format(mix_vle.cp_mass))