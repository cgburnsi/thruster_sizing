from unbar import UNBAR

# Initial Concentration of Species in the Vapor Region of the Reactor Bed
def CONC(T, P, WM4, WM3, WM2, WM1, R, H, HF):
    TVAP, PRES = T, P
    H4 = UNBAR('TBLH4', TVAP)
    XV = -(H - HF) / H4
    C4 = (PRES * WM4) / (R * TVAP) * ((1.0 - XV) / (1.0 + XV))          # Hydrazine Concentration
    C3 = (PRES * WM3) / (R * TVAP) * (XV / (1.0 + XV))                  # Ammonia Concentration
    C2 = (PRES * WM2) / (2.0 * R * TVAP) * (XV / (1.0 + XV))            # Nitrogen Concentration
    C1 = (PRES * WM1) / (2.0 * R * TVAP) * (XV / (1.0 + XV))            # Hydrogen Concentration
    return C1, C2, C3, C4
