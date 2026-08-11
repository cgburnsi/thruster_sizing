
class TemperatureInterval:
    """Holds data for a single temperature range (Records 3, 4, 5)."""
    def __init__(self, temp_low, temp_high, cp_coeffs, b1, b2):
        self.temp_low = temp_low
        self.temp_high = temp_high
        self.cp_coeffs = cp_coeffs
        self.b1 = b1
        self.b2 = b2
        
    def __repr__(self):
        """Gives a nice print-out for debugging."""
        return f"Interval({self.temp_low}K - {self.temp_high}K)"
        

class Species:
    """Holds all data for one chemical species."""
    def __init__(self, name, comments, phase, molecular_weight, hf298, chemical_formula):
        # Data from Records 1 & 2
        self.name = name
        self.comments = comments
        self.phase = phase
        self.molecular_weight = molecular_weight
        self.hf298 = hf298
        self.chemical_formula = chemical_formula
        
        # We initialize an empty list to hold the interval objects
        self.intervals = []

    def add_interval(self, interval_object):
        """A helper method to add a new interval."""
        self.intervals.append(interval_object)

    def __repr__(self):
        """Gives a nice print-out for debugging."""
        return f"Species(name='{self.name}', phase={self.phase}, intervals={len(self.intervals)})"
    
    


# --- This is the full example of parsing one species ---

if __name__ == '__main__':

    # This dictionary will hold all your species objects
    all_species_db = {}
    
    # --- Example data for one species (ALBO2) ---
    # In your real code, you would be reading these lines from a file
    line1 = "ALBO2             JANAF 6/66 JPCRD v14 sup 1 1985.                              "
    line2 = ' 2 J 6/66 AL  1.00B   1.00O   2.00    0.00    0.00 0   69.7913380    -541416.208'
    # First temperature interval
    line3 = '    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000'
    line4 = ' 0.000000000D+00 0.000000000D+00 2.308723400D+00 1.889053900D-02-2.063334800D-05'
    line5 = ' 1.025132400D-08-1.694128300D-12                -6.648216700D+04 1.447701848D+01'
    # Second temperature interval
    line6 = '   1000.000   5000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000'
    line7 = '-1.579040000D+06 8.361595000D+03-2.012111000D+00 3.197931000D-03-3.059902000D-06'
    line8 = ' 8.079231000D-10-8.211904000D-14                -6.559810100D+04-1.077202300D+01'
    # (Note: I've added a dummy set of lines 6-8 for the 2nd interval)
    
    # --- PARSE RECORD 1 ---
    species_name = line1[0:16].strip()
    comments = line1[18:80].strip()

    # --- PARSE RECORD 2 ---
    num_T_intervals = int(line2[0:2].strip())
    phase = int(line2[50:52].strip())
    mw = float(line2[52:65].strip())
    hf298 = float(line2[65:80].strip())

    # Parse chemical formula
    formula = []
    for i in range(5):
        start_atom = 10 + (i * 8)
        start_num = 12 + (i * 8)
        atom = line2[start_atom : start_num].strip()
        if atom:
            num = float(line2[start_num : start_num + 6].strip())
            formula.append((atom, num))

    # --- CREATE THE SPECIES OBJECT ---
    # Now we can create the main object
    current_species = Species(
        name=species_name,
        comments=comments,
        phase=phase,
        molecular_weight=mw,
        hf298=hf298,
        chemical_formula=formula
    )

    # --- LOOP FOR EACH TEMPERATURE INTERVAL ---
    # Now, your code must loop 'num_T_intervals' times
    # This example is hard-coded for 2 intervals
    
    # --- Process Interval 1 (Lines 3, 4, 5) ---
    t_low1 = float(line3[0:11].strip())
    t_high1 = float(line3[11:22].strip())
    
    # (Don't forget to replace 'D' with 'E' for Python float conversion)
    coeffs1 = [
        float(line4[0:16].strip().replace('D', 'E')),
        float(line4[16:32].strip().replace('D', 'E')),
        float(line4[32:48].strip().replace('D', 'E')),
        float(line4[48:64].strip().replace('D', 'E')),
        float(line4[64:80].strip().replace('D', 'E')),
        float(line5[0:16].strip().replace('D', 'E')),
        float(line5[16:32].strip().replace('D', 'E'))
    ]
    b1_const1 = float(line5[48:64].strip().replace('D', 'E'))
    b2_const1 = float(line5[64:80].strip().replace('D', 'E'))
    
    # Create the interval object
    interval1 = TemperatureInterval(t_low1, t_high1, coeffs1, b1_const1, b2_const1)
    
    # Add it to the species
    current_species.add_interval(interval1)


    # --- Process Interval 2 (Lines 6, 7, 8) ---
    t_low2 = float(line6[0:11].strip())
    t_high2 = float(line6[11:22].strip())
    
    coeffs2 = [
        float(line7[0:16].strip().replace('D', 'E')),
        float(line7[16:32].strip().replace('D', 'E')),
        float(line7[32:48].strip().replace('D', 'E')),
        float(line7[48:64].strip().replace('D', 'E')),
        float(line7[64:80].strip().replace('D', 'E')),
        float(line8[0:16].strip().replace('D', 'E')),
        float(line8[16:32].strip().replace('D', 'E'))
    ]
    b1_const2 = float(line8[48:64].strip().replace('D', 'E'))
    b2_const2 = float(line8[64:80].strip().replace('D', 'E'))
    
    interval2 = TemperatureInterval(t_low2, t_high2, coeffs2, b1_const2, b2_const2)
    current_species.add_interval(interval2)

    
    # --- FINAL STEP ---
    # Add the fully-parsed species to your main database
    all_species_db[current_species.name] = current_species

    # --- NOW YOU CAN EASILY ACCESS THE DATA ---
    print(f"Successfully parsed: {all_species_db['ALBO2']}")
    print(f"Molecular weight: {all_species_db['ALBO2'].molecular_weight}")
    print(f"Chemical formula: {all_species_db['ALBO2'].chemical_formula}")
    print("\nIntervals:")
    print(f"  Interval 1: {all_species_db['ALBO2'].intervals[0]}")
    print(f"  Interval 2: {all_species_db['ALBO2'].intervals[1]}")
    print(f"\nCoeff 'a3' from Interval 1: {all_species_db['ALBO2'].intervals[0].cp_coeffs[2]}")
    print(f"Coeff 'b1' from Interval 2: {all_species_db['ALBO2'].intervals[1].b1}")
    
    
    
    
    
    
    
    
    