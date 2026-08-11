
from thermoprop import SpeciesDB



import math

R_UNIV = 8.314462618  # J/mol-K

def _pick_range(sp_record, T):
    for r in sp_record['thermo_ranges']:
        if r['t_low'] <= T <= r['t_high']:
            return r
    raise ValueError(f'no NASA7 range covers T={T} K for {sp_record.get("species","?")}')

def nasa7_cp_h_s(sp_record, T):
    r = _pick_range(sp_record, T)
    # support either 'a' (10 terms) or 'coefficients' (7 terms) if you ever used that name
    coeffs = r.get('a') or r.get('nasa7')
    if coeffs is None:
        raise KeyError('range missing coefficients')
    # If you ever store only 7 coeffs, pad to 10 for uniform handling
    a = list(coeffs) + [0.0] * (10 - len(coeffs))

    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = a
    t = T

    # J/mol-K
    cp_molar = R_UNIV*(a1 + a2*t + a3*t**2 + a4*t**3 + a5*t**4)

    # J/mol
    h_molar  = R_UNIV*t*(a1 + a2*t/2 + a3*t**2/3 + a4*t**3/4 + a5*t**4/5 + a6/t)

    # J/mol-K
    s_molar  = R_UNIV*(a1*math.log(t) + a2*t + a3*t**2/2 + a4*t**3/3 + a5*t**4/4 + a7)

    return cp_molar, h_molar, s_molar

def cp_from_json_record(sp_record, T, basis='molar'):
    cp_molar, _, _ = nasa7_cp_h_s(sp_record, T)
    if basis.lower() == 'molar':
        return cp_molar  # J/mol-K
    # mass basis
    mw_g_per_mol = sp_record.get('mw')
    if mw_g_per_mol is None:
        # fallback if you used a different key in some files
        mw_g_per_mol = sp_record.get('molar_mass_g_mol')
    if mw_g_per_mol is None:
        raise KeyError('species record missing molar mass (mw)')
    mw = mw_g_per_mol / 1000.0  # kg/mol
    return cp_molar / mw  # J/kg-K




if __name__ == '__main__':
    
    db = SpeciesDB(validate=True)
    
    db.list_species()

    n2 = db.get('N2')

    T = 300  # K
    cp_molar = cp_from_json_record(n2, T, basis='molar')
    cp_mass  = cp_from_json_record(n2, T, basis='mass')
    
    print('Cp_molar(N2, 1200 K) =', cp_molar, 'J/mol-K')
    print('Cp_mass (N2, 1200 K) =', cp_mass,  'J/kg-K')
        