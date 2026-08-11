def print_inputs(HF, HL, HV, TF, TVAP, CFL, PRESSURE, KP, F, G0, R, ALPHA3, CGM, DIF3, DIF4, WM4, WM3, WM2, WM1, ZEND, AGM, BGM, ALPHA1, ALPHA2, N1, N2, N3, ENMX1, ENMX2, ENMX3, Z0):
    header_list = ['HF', 'HL', 'HV', 'TF', 'TVAP', 'CFL', 'PRESSURE', 'KP', 'F', 'G0', 'R', 'ALPHA3', 'CGM', 'DIF3', 'DIF4', 'WM4', 'WM3', 'WM2', 'WM1', 'ZEND', 'AGM', 'BGM', 'ALPHA1', 'ALPHA2', 'N1', 'N2', 'N3', 'ENMX1', 'ENMX2', 'ENMX3', 'Z0']
    # find max length
    head_len = [len(i) for i in header_list]
    max_head = max(head_len)

    # pad blanks
    pretty_head = [i.ljust(max_head) for i in header_list]
    value_list = [HF, HL, HV, TF, TVAP, CFL, PRESSURE, KP, F, G0, R, ALPHA3, CGM, DIF3, DIF4, WM4, WM3, WM2, WM1, ZEND, AGM, BGM, ALPHA1, ALPHA2, N1, N2, N3, ENMX1, ENMX2, ENMX3, Z0]
    value_strings = [str(i) for i in value_list]
    value_len = [len(i) for i in value_strings]
    max_value = max(value_len)

    # pad blanks
    pretty_value = [i.ljust(max_value) for i in value_strings]
    
    unit_list = ['[BTU/lb]', '[BTU/lb]', '[BTU/lb]', '[degR]', '[degR]', '[degR]]', '[psia]', '[[Btu/ft-sec-degR]', '[lb/ft3-sec]', 
                 '[lb/ft2-s]', '[psia-ft3/lb-mol-degR]', '[1/sec]', '[degR]', '[ft2/s]', '[ft2/s]', '[lb/lb-mol]', '[lb/lb-mol]', 
                 '[lb/lb-mol]', '[lb/lb-mol]', '[ft]', '[degR]', '[degR]', '[1/sec]', '[(lb/ft^3)^1.6/sec]', '[-]', '[-]', '[-]', '[-]', 
                 '[-]', '[-]', '[ft]']
    # find max length
    unit_len = [len(i) for i in unit_list]
    max_unit = max(unit_len)

    # pad blanks
    pretty_unit = [i.ljust(max_unit) for i in unit_list]
    
    # merge strings
    pretty_front = [" = ".join(i) for i in zip(pretty_head, pretty_value)]
    pretty_pretty = [" ".join(i) for i in zip(pretty_front, pretty_unit)]
    pretty_pretty.insert(0, "=" * (max_head + max_value + max_unit + 4))
    pretty_pretty.append("=" * (max_head + max_value + max_unit + 4))
    
    print('\nInput Constants')
    for p in pretty_pretty: print(p)

def print_liquid_data(Z, TEMP, H, DHDZ):
    print(f'Z:{Z:.3e} [ft], TEMP:{TEMP:.2f} [degR], H:{H:.2f} [BTU/lb], DHDZ:{DHDZ:.3e} [ft]')
    
def print_lqvp_data(Z, TEMP, H, WFV):
    print(f'Z:{Z:.3e} [ft], TEMP:{TEMP:.2f} [degR], H:{H:.2f} [BTU/lb], WFV:{WFV:14.5E}')

# --- NEW FUNCTION ---
def print_vapor_data(Z, TEMP, PRES, H, C1, C2, C3, C4, FRAC3D):
    print(f'Z:{Z:.6e} [ft], T:{TEMP:.2f} [R], P:{PRES:.2f} [psia], H:{H:.2f} [BTU/lb], '
          f'C4:{C4:.3e}, C3:{C3:.3e}, C2:{C2:.3e}, C1:{C1:.3e}, FRAC3D:{FRAC3D:.4f}')

# --- NEW DEBUG FUNCTION ---
def print_vapor_debug(Z, T, DHDZ, DZ, T2, T3, T4, H3, H4, TGRAD, HC_term, F_term):
    print(f"  Z={Z:.4e} | T={T:.2f} | DHDZ={DHDZ:.3e} | DZ={DZ:.3e}")
    print(f"    Rates (lb/ft^3-s):")
    print(f"      T2 (N2H4_therm) = {T2:.3e}")
    print(f"      T3 (NH3_cat)    = {T3:.3e}")
    print(f"      T4 (N2H4_cat)   = {T4:.3e}")
    print(f"    Heat Terms (Btu/ft^3-s):")
    print(f"      N2H4_Heat = -H4*(T2+T4) = {-H4*(T2+T4):.3e}")
    print(f"      NH3_Heat  = -H3*T3      = {-H3*T3:.3e}")
    print(f"      HC_Loss   = -HC*AP*(T-TGRAD) = {HC_term:.3e}")
    print(f"      F_Loss    = -F*(H-HF)   = {F_term:.3e}")
    print(f"    (H4={H4:.1f}, H3={H3:.1f}, TGRAD={TGRAD:.2f})")