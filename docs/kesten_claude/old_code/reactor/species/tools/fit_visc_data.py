
import numpy as np

def fit_mu_coeffs(T, mu_Pa_s, form='mu_log10_microP_logT', weights=None):
    """
    Fit viscosity correlation coefficients from data.
    T: array-like [K]
    mu_Pa_s: array-like [Pa·s]
    form: 'mu_log10_microP_logT' (4-term) or 'mu_log10_microP' (3-term)
    weights: optional array-like; higher = more influence on those points
    Returns tuple of coefficients (A,B,C,D) or (A,B,C)
    """
    T = np.asarray(T, float)
    mu = np.asarray(mu_Pa_s, float)
    if np.any(T <= 0) or np.any(mu <= 0):
        raise ValueError("All T and μ must be > 0.")

    # target: y = log10(mu in microPoise) = log10(mu / 1e-7)
    y = np.log10(mu / 1e-7)

    if form == 'mu_log10_microP_logT':
        X = np.column_stack([np.ones_like(T), 1.0/T, 1.0/(T*T), np.log10(T)])
    elif form == 'mu_log10_microP':
        X = np.column_stack([np.ones_like(T), 1.0/T, 1.0/(T*T)])
    else:
        raise ValueError("Unknown form: %r" % form)

    if weights is not None:
        w = np.sqrt(np.asarray(weights, float))
        Xw = X * w[:, None]
        yw = y * w
        coeffs, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    else:
        coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)

    return tuple(coeffs.tolist())


def mu_from_coeffs(T, coeffs, form='mu_log10_microP_logT'):
    T = np.asarray(T, float)
    if form == 'mu_log10_microP_logT':
        A,B,C,D = coeffs
        y = A + B/T + C/(T*T) + D*np.log10(T)
    elif form == 'mu_log10_microP':
        A,B,C = coeffs
        y = A + B/T + C/(T*T)
    else:
        raise ValueError("Unknown form.")
    # convert back: μ[Pa·s] = 10^y × 1e-7
    return (10.0**y) * 1e-7


if __name__ == '__main__':
    
    # Your experimental/simulated data
    T_data  = np.array([180.0,   360.0,  534.6,   540.0,   720.0,   900.0,  1080.0], float)
    mu_data = np.array([1390.16, 1332.82, 1280.2, 1279.12, 1237.79, 1208.80, 1189.76], float)  # Pa·s (example)
    
    # Fit
    A,B,C,D = fit_mu_coeffs(T_data, mu_data, form='mu_log10_microP_logT')
    print("A,B,C,D =", A,B,C,D)
    
    # Check fit quality
    mu_fit = mu_from_coeffs(T_data, (A,B,C,D), form='mu_log10_microP_logT')
    rel_err = (mu_fit - mu_data)/mu_data
    print("Max |rel err|:", np.max(np.abs(rel_err)))



'''
   elif table.upper() == 'DHVST':
       T     =   [180.0,   360.0,  534.6,   540.0,   720.0,   900.0,  1080.0]
       DHVST = [1390.16, 1332.82, 1280.2, 1279.12, 1237.79, 1208.80, 1189.76]
       calc_val = np.interp(data_in, T, DHVST)  # [?] 

   elif table.upper() == 'DHLVST':
       T      = [ 180.0,  360.0, 534.6,   540.0,  720.0,  900.0, 1080.0]
       DHLVST = [652.14, 665.96, 679.61, 679.89, 700.89, 733.19, 777.22]
       calc_val = np.interp(data_in, T, DHLVST)  # [?]

   else:
       print(f'ERROR - table: {table}, value: {data_in} not calculated correctly')
       calc_val = -1
   return calc_val

'''


