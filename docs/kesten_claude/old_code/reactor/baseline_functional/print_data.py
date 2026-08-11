
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
