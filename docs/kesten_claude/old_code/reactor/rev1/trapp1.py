import numpy as np

def RHETF(A, B, C, D, E, N):
    return E*A**(1-N)*B**N * np.exp(C*D*(1.0 - B/A)/(1+D*(1-B/A)))    

def FOXI1(X, R):
    return X**2 * R
    
def CPXF(X, Y, Z):
    return (X-Y)/(1-Y) * Z
 
# Trapelzoidal Method Integration
def TRAPP(U, V, NPART, XOA, CPS, GAMMA, BETA, K0):
    N = NPART - 1
    PART = NPART
    H = (U - V) / PART
    UPH = U + H
    SUM = 0.0
    CPX1 = CPXF(U, XOA, CPS)
    CPX2 = CPXF(V, XOA, CPS)
    RHET1 = RHETF(CPS, CPX1, GAMMA, BETA, K0, 1)
    RHET2 = RHETF(CPS, CPX2, GAMMA, BETA, K0, 1)
    
    # Calculate first and last terms of riemann sum first
    TRM1 = FOXI1(U, RHET1)/2.0
    TRM2 = FOXI1(V, RHET2)/2.0
    
    # Do the rule
    for idx in range(1,N):
        CPX = CPXF(UPH, XOA, CPS)
        RHET = RHETF(CPS, CPX, GAMMA, BETA, K0, 1)
        SUM = SUM + FOXI1(UPH, RHET)
        UPH = UPH + H
        
    RIESUM = H * (TRM1 + SUM + TRM2)

    return RIESUM

