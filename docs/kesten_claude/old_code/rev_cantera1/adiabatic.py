


import cantera as ct
import numpy as np



if __name__ == '__main__':
    
    g = ct.Solution('rxn_mechanisms.yaml', 'gas')
    
    # Example: 2 parts H2, 1 part H2O, and 4 parts He
    g.TPX = 1000.0, ct.one_atm, 'H2:2.0, H2O:1.0, He:4.0'

    g.equilibrate('HP')
    
    rf = g.forward_rates_of_progress
    rr = g.reverse_rates_of_progress
    for i in range(g.n_reactions):
        if g.reaction(i).reversible and rf[i] != 0.0:
            print(' %4i  %10.4g  ' % (i, (rf[i] - rr[i])/rf[i]))



