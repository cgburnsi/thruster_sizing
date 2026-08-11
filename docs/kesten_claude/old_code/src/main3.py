
class TempInterval:
    def __init__(self, T_low, T_high, coeffs, b1, b2):
        self.T_low = T_low
        self.T_high = T_high
        self.coeffs = coeffs
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
    
    


if __name__ == '__main__':
        
    # This is the example line for Record 1 from the document
    # It's important that this string is *exactly* 80 characters long
    # and matches the file format.
    
    line = "ALBO2             JANAF 6/66 JPCRD v14 sup 1 1985.                              "
    # --- New, More Robust Method ---
    
    # 1. Parse the Species from its *defined* columns (1-16)
    #    Slice: line[0:16]
    species_name = line[0:16].strip()
    
    # 2. Parse the Comments. We'll *still* use the spec (19-80),
    #    but we must be sure the line data is correct.
    #    Slice: line[18:80]
    #    Column 19 is index 18.
    comments = line[18:80].strip()
    
    
    line2 = ' 2 J 6/66 AL  1.00B   1.00O   2.00    0.00    0.00 0   69.7913380    -541416.208'
    num_T_intervals = line2[0:2].strip()
    ref_date_code   = line2[3:9].strip()
    species1        = line2[10:12].strip()
    s1num           = line2[12:18].strip()
    species2        = line2[18:20].strip()
    s2num           = line2[20:26].strip()
    species3        = line2[26:28].strip()
    s3num           = line2[28:34].strip()
    species4        = line2[34:36].strip()
    s4num           = line2[36:42].strip()
    species5        = line2[42:44].strip()
    s5num           = line2[44:50].strip()
    phase           = line2[50:52].strip()
    mw              = line2[52:65].strip()
    hf298           = line2[65:80].strip()
    
    line3 = '    300.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000'
    temp_range_low  = line3[0:11]
    temp_range_high = line3[11:22]
    num_ceoff       = line3[22:23]
    temp_exponents1 = line3[23:28]
    temp_exponents2 = line3[28:33]
    temp_exponents3 = line3[33:38]
    temp_exponents4 = line3[38:43]
    temp_exponents5 = line3[43:48]
    temp_exponents6 = line3[48:53]
    temp_exponents7 = line3[53:58]
    temp_exponents8 = line3[58:63]
    deltaH298       = line3[65:80]
    
    line4 = ' 0.000000000D+00 0.000000000D+00 2.308723400D+00 1.889053900D-02-2.063334800D-05'
    cp_coeff1       = line4[0:16]
    cp_coeff2       = line4[16:32]
    cp_coeff3       = line4[32:48]
    cp_coeff4       = line4[48:64]
    cp_coeff5       = line4[64:80]
    
    line5 = ' 1.025132400D-08-1.694128300D-12                -6.648216700D+04 1.447701848D+01'
    cp_coeff6       = line5[0:16]
    cp_coeff7       = line5[16:32]
    b1              = line5[48:64]
    b2              = line5[64:80]
    
    
    
    
    
    
    
    