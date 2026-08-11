
class TempRange:
    def __init__(self, T_min=0.0, T_max=0.0, coeffs=None):
        self.T_min  = T_min
        self.T_max  = T_max
        
        if coeffs is None: 
            self.coeffs = [] 
        else: 
            self.coeffs = coeffs
            
    def __repr__(self):
        return f"TempRange(t_min={self.t_min}, t_max={self.t_max}, coeffs=[...{len(self.coeffs)} coeffs...])"

class ThermoEntry:
    def __init__(self, name="", source_info="", elements=None, mw=0.0, hf298=0.0, ranges=None):
        self.name = name
        self.source_info = source_info
        if elements is None:
            self.elements = {}
        else:
            self.elements = elements
        self.mw = mw
        self.hf298 = hf298
        if ranges is None:
            self.ranges = []
        else:
            self.ranges = ranges
            
    # I've added a __repr__ so it still prints nicely like a dataclass
    def __repr__(self):
        parts = []
        parts.append(f"name='{self.name}'")
        parts.append(f"source_info='{self.source_info}'")
        parts.append(f"elements={self.elements}")
        parts.append(f"mw={self.mw}")
        parts.append(f"hf298={self.hf298}")
        parts.append(f"ranges={self.ranges}")
        return f"ThermoEntry({', '.join(parts)})"

if __name__ == '__main__':
    
    # Variables from the previous example
    STATE_PRE_THERMO = 0
    STATE_PRODUCTS = 1
    STATE_REACTANTS = 2
    
    current_state = STATE_PRE_THERMO
    current_record_key = None
    data = {"products": {}, "reactants": {}}
    
    with open('nasa9.dat') as f:
        
        # Variables from the previous example
        STATE_PRE_THERMO = 0
        STATE_PRODUCTS = 1
        STATE_REACTANTS = 2
        
        current_state = STATE_PRE_THERMO
        current_record_key = None
        data = {"products": {}, "reactants": {}}
        
        # --- Inside your 'with open(...)' loop ---
        for line in f:
            processed_line = line.strip()
            
            # 1. Skip comments/blank lines
            if not processed_line or processed_line.startswith('!'):
                continue
        
            # 2. Check for state-changing keywords
            if processed_line == 'thermo':
                current_state = STATE_PRODUCTS
                current_record_key = None # Reset the record key
                continue # Move to the next line
        
            elif processed_line == 'END PRODUCTS':
                current_state = STATE_REACTANTS
                current_record_key = None # Reset the record key
                continue # Move to the next line
        
            elif processed_line == 'END REACTANTS':
                break # We are done
        
            # 3. If it's not a keyword or comment, process it based on the state
            # This is the part your code was missing.
            
            if current_state == STATE_PRODUCTS:
                if processed_line[0].isalpha():
                    current_record_key = processed_line
                    data["products"][current_record_key] = []
                elif current_record_key:
                    # This is numerical data for the current product
                    data["products"][current_record_key].append(processed_line.split())
        
            elif current_state == STATE_REACTANTS:
                if processed_line[0].isalpha():
                    current_record_key = processed_line
                    data["reactants"][current_record_key] = []
                elif current_record_key:
                    # This is numerical data for the current reactant
                    data["reactants"][current_record_key].append(processed_line.split())
                    
            # Note: If current_state is STATE_PRE_THERMO,
            # any data lines here are just ignored, which is probably what you want.
            
