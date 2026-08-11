
class TemperatureInterval:
    """Holds data for a single temperature range (Records 3, 4, 5)."""
    def __init__(self, temp_low, temp_high, cp_expons, cp_coeffs, b1, b2):
        self.temp_low = temp_low
        self.temp_high = temp_high
        self.cp_expons = cp_expons
        self.cp_coeffs = cp_coeffs
        self.b1 = b1
        self.b2 = b2
        
    def __repr__(self):
        """Gives a nice print-out for debugging."""
        return f"Interval({self.temp_low}K - {self.temp_high}K)"

class Species:
    """Holds all data for one chemical species."""
    def __init__(self, name, comments, phase, molecular_weight, hf298, chemical_formula, role):
        self.name = name
        self.comments = comments
        self.phase = phase
        self.molecular_weight = molecular_weight
        self.hf298 = hf298
        self.chemical_formula = chemical_formula
        self.role = role
        self.intervals = []

    def add_interval(self, interval_object):
        """A helper method to add a new interval."""
        self.intervals.append(interval_object)

    def __repr__(self):
        """Gives a nice print-out for debugging."""
        return f"Species(name='{self.name}', role='{self.role}', intervals={len(self.intervals)})"

class Nasa9Parser:
    """
    Parses NASA-9 format thermodynamic data files by bundling all
    schemas and helper logic into a single class.
    
    Usage:
        parser = Nasa9Parser()
        db = parser.parse('nasa9.dat')
    """

    # --- Static Helper Functions ---
    # These are marked '@staticmethod' because they don't need 'self'.
    # They are pure utility functions.
    @staticmethod
    def _parse_str(val):     return val.strip()
    
    @staticmethod
    def _parse_int(val):     return int(val.strip())
    
    @staticmethod
    def _parse_D_float(val): return float(val.replace('D', 'E'))
    
    @staticmethod
    def _parse_float(val):   return 0.0 if not val.strip() else float(val.strip())

    @staticmethod
    def _assemble_formula(record_data):
        formula = []
        for i in range(1, 6):
            atom_key, num_key = f'element{i}_name', f'element{i}_num'
            atom = record_data[atom_key]
            if atom: formula.append((atom, record_data[num_key]))
        return formula
    
    @staticmethod
    def _parse_line(line, schema):
        data = {}
        for name, (line_slice, parse_function) in schema.items():
            # 'parse_function' refers to one of the static methods above
            data[name] = parse_function(line[line_slice])
        return data

    # --- NASA9 Data File Parsing Schemas ---
    # These are class variables. They refer to the static methods defined above.
    RECORD_1_SCHEMA = {'name':          (slice( 0, 16), _parse_str),
                       'comments':      (slice(18, 80), _parse_str)}
    RECORD_2_SCHEMA = {'num_intervals': (slice( 0,  2), _parse_int),
                       'element1_name': (slice(10, 12), _parse_str), 'element1_num': (slice(12, 18), _parse_float),
                       'element2_name': (slice(18, 20), _parse_str), 'element2_num': (slice(20, 26), _parse_float),
                       'element3_name': (slice(26, 28), _parse_str), 'element3_num': (slice(28, 34), _parse_float),
                       'element4_name': (slice(34, 36), _parse_str), 'element4_num': (slice(36, 42), _parse_float),
                       'element5_name': (slice(42, 44), _parse_str), 'element5_num': (slice(44, 50), _parse_float),                  
                       'phase':         (slice(50, 52), _parse_int),
                       'mw':            (slice(52, 65), _parse_float),
                       'hf298':         (slice(65, 80), _parse_float)}
    RECORD_3_SCHEMA = {'t_low':         (slice( 0, 11),  _parse_float), 't_high': (slice(11, 22), _parse_float),
                       'num_ceoff':     (slice(22, 23), _parse_int),
                       'e1':            (slice(23, 28), _parse_float),
                       'e2':            (slice(28, 33), _parse_float),
                       'e3':            (slice(33, 38), _parse_float),
                       'e4':            (slice(38, 43), _parse_float),
                       'e5':            (slice(43, 48), _parse_float),
                       'e6':            (slice(48, 53), _parse_float),
                       'e7':            (slice(53, 58), _parse_float),
                       'e8':            (slice(58, 63), _parse_float),
                       'deltaH298':     (slice(65, 80), _parse_float)}
    RECORD_4_SCHEMA = {'a1': (slice( 0, 16), _parse_D_float), 'a2': (slice(16, 32), _parse_D_float),
                       'a3': (slice(32, 48), _parse_D_float), 'a4': (slice(48, 64), _parse_D_float),
                       'a5': (slice(64, 80), _parse_D_float)}
    RECORD_5_SCHEMA = {'a6': (slice( 0, 16), _parse_D_float), 'a7': (slice(16, 32), _parse_D_float),
                       'b1': (slice(48, 64), _parse_D_float), 'b2': (slice(64, 80), _parse_D_float)}


    @staticmethod
    def _parse_and_add_intervals(current_species, f_iter, num_intervals):
        """
        Parses and adds all temperature intervals for a given species.
        Also handles the '0-interval' case by skipping the extra line.
        """
        if num_intervals > 0:
            for _ in range(num_intervals):
                record_3 = Nasa9Parser._parse_line(next(f_iter), Nasa9Parser.RECORD_3_SCHEMA)
                record_4 = Nasa9Parser._parse_line(next(f_iter), Nasa9Parser.RECORD_4_SCHEMA)
                record_5 = Nasa9Parser._parse_line(next(f_iter), Nasa9Parser.RECORD_5_SCHEMA)       
        
                expons = [record_3['e1'], record_3['e2'], record_3['e3'], record_3['e4'],
                          record_3['e5'], record_3['e6'], record_3['e7'], record_3['e8']]
                coeffs = [record_4['a1'], record_4['a2'], record_4['a3'],record_4['a4'],
                          record_4['a5'],record_5['a6'], record_5['a7']]
            
                interval = TemperatureInterval(temp_low  = record_3['t_low'],
                                               temp_high = record_3['t_high'],
                                               cp_expons = expons,
                                               cp_coeffs = coeffs,
                                               b1 = record_5['b1'], b2 = record_5['b2'])
                
                current_species.add_interval(interval)
        else:
            # This is the 0-interval case (like B2H6(L))
            # We must read and discard its single 'Record 3' line
            next(f_iter)

    # --- Main Public Method ---
    # This is the single function you call to run the parser.
    def parse(self, datafile):
        """
        Parses a NASA-9 format thermodynamic data file.
        
        Args:
            datafile (str): The path to the data file (e.g., "nasa9.dat")
            
        Returns:
            dict: A dictionary of {species_name: Species_Object}
        """
        species_db  = {}
        current_role = 'unknown'
        
        with open(datafile, 'r') as f:
            f_iter = iter(f)
    
            while True:
                try: line1 = next(f_iter)
                except StopIteration: break
                
                # 1. Handle all header and comment lines
                line_clean = line1.strip()
                if line_clean.startswith("thermo"):
                    current_role = "product"                
                    next(f_iter) # Skip global temp line
                    continue 
                if line_clean.startswith("END PRODUCTS"):
                    current_role = "reactant"
                    continue
                if line_clean.startswith("END"): break
                if line_clean.startswith("!"): continue
            
                # 2. Parse Species Records 1 & 2
                # We call the static methods on the class itself
                record_1 = Nasa9Parser._parse_line(line1, Nasa9Parser.RECORD_1_SCHEMA)
                record_2 = Nasa9Parser._parse_line(next(f_iter), Nasa9Parser.RECORD_2_SCHEMA)     
                
                # 3. Create the Species object
                formula = Nasa9Parser._assemble_formula(record_2)
                current_species = Species(name  = record_1['name'],  comments = record_1['comments'],
                                          phase = record_2['phase'], molecular_weight = record_2['mw'],
                                          hf298 = record_2['hf298'], chemical_formula=formula,
                                          role  = current_role)
                
                # 4. Parse and Add all its intervals
                Nasa9Parser._parse_and_add_intervals(current_species, f_iter, record_2['num_intervals'])
                    
                # 5. Add the finished species to the database
                species_db[current_species.name] = current_species
            
        return species_db


# --- This is the main part that runs your code ---
if __name__ == '__main__':
    
    print("Parsing nasa9.dat...")
    
    # 1. Create an instance of the parser
    parser = Nasa9Parser()
    
    # 2. Call the 'parse' method on that instance
    db = parser.parse('nasa9.dat')
    
    print(f"--- Parsing Complete ---")
    print(f"Successfully parsed {len(db)} species.")

    # 3. Print some examples
    if "H2O" in db:
        print(f"\nExample (product): {db['H2O']}")
        print(f"  Formula: {db['H2O'].chemical_formula}")
        print(f"  Num Intervals: {len(db['H2O'].intervals)}")

    if "B2H6(L)" in db:
        print(f"\nExample (reactant): {db['B2H6(L)']}")
        print(f"  Formula: {db['B2H6(L)'].chemical_formula}")
        print(f"  Num Intervals: {len(db['B2H6(L)'].intervals)}")