import json
import numpy as np
from pathlib import Path

# ---------- FLUID PROPERTIES -------------------------------------------------------------------------------
class PropertyTables:

    def __init__(self, json_path='data.json'):
        self.tables = {}
        self.path   = Path(json_path)
        
        if not self.path.exists(): raise FileNotFoundError(f'Could not find properties file: {json_path}')

        with open(self.path, 'r') as f:
            data = json.load(f)
            
        for name, table_data in data['tables'].items():
            self.tables[name] = self._parse_table(table_data)
        
        # PATCH: Restore high-precision/legacy data from the original Kesten code
        # where the JSON file was truncated or inaccurate.
        self._apply_kesten_patches()
 
            
    def _parse_table(self, table_dict):
        keys  = list(table_dict.keys())
        x_key = next(k for k in keys if k not in ["y", "description"])
        
        return np.array(table_dict[x_key]), np.array(table_dict["y"])
    
    def _interp(self, table_name, value):
        if table_name not in self.tables:
            raise KeyError(f'Table "{table_name}" Not Found In Loaded Fluid Properties')
            
        x, y = self.tables[table_name]
        return np.interp(value, x, y)
    
    def _apply_kesten_patches(self):
        """
        Hardcoded corrections to restore the specific datasets used in the original 
        UNBAR routine. TBLVP and TVAP (VPTBL) are distinct datasets.
        """
        # 1. TBLVP: Temperature [R] -> Vapor Pressure [psia]
        # Used for local vapor pressure in the liquid region.
        t_vp = [492.0, 519.0, 528.37, 529.08, 534.60, 534.71, 538.84, 543.91, 
                545.73, 560.20, 569.98, 579.26, 579.48, 595.34, 610.13, 614.08, 
                618.07, 627.49, 628.82, 645.68, 650.76, 665.57, 674.99 ,686.13, 
                692.39, 697.47, 744.0, 798.0, 852.0, 942.0, 1032.0, 1122.0, 
                1176.0]
        p_vp = [0.0520, 0.1479, 0.2011, 0.2069, 0.2398, 0.2436, 0.2823, 
                0.2920, 0.3539, 0.5453, 0.7367, 0.9727, 0.9823, 1.510, 
                2.204, 2.462, 2.740, 3.407, 3.562, 5.240, 5.971, 8.065, 
                9.711, 11.91, 13.46, 14.70, 33.80, 73.48, 147.0, 382.1, 
                823.0, 1528.0, 2131.0]

        self.tables['TBLVP'] = (np.array(t_vp), np.array(p_vp))

        # 2. VPTBL (aka TVAP): Pressure [psia] -> Saturation Temp [R]
        # Used for finding Tsat from System Pressure. 
        # NOTE: This uses a DIFFERENT set of points than TBLVP in the original code.
        p_tvap = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0,
                 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0,
                 850.0, 900.0, 950.0, 1000.0]
        t_tvap = [770.0, 820.0, 855.0, 880.0, 905.0, 925.0, 945.0, 965.0,
                980.0, 995.0, 1010.0, 1025.0, 1035.0, 1050.0, 1060.0, 
                1070.0, 1080.0, 1090.0, 1100.0, 1110.0]
        
        # Ensure x=Pressure, y=Temperature
        self.tables['VPTBL'] = (np.array(p_tvap), np.array(t_tvap))

        # 3. VISVST: Viscosity correction
        t_vis = [360.0, 540.0, 720.0, 900.0, 1080.0, 1260.0, 1440.0, 1620.0, 
                 1800.0, 1980.0, 2160.0, 2340.0, 2520.0]
        y_vis = [0.048e-4, 0.070e-4, 0.093e-4, 0.117e-4, 0.141e-4, 0.164e-4, 
                 0.186e-4, 0.207e-4, 0.220e-4, 0.247e-4, 0.266e-4, 0.285e-4, 
                 0.302e-4]
        self.tables['VISVST'] = (np.array(t_vis), np.array(y_vis))

    # --- Methods ---
    def H3TBL(self, X):  return self._interp('H3TBL', X)
    def H4TBL(self, X):  return self._interp('H4TBL', X)
    def DHVST(self, X):  return self._interp('DHVST', X)
    def DHLVST(self, X): return self._interp('DHLVST', X)
    def VISVST(self, X): return self._interp('VISVST', X)
    def CFTBL1(self, X): return self._interp('CFTBL1', X)
    def CFTBL2(self, X): return self._interp('CFTBL2', X)
    def CFTBL3(self, X): return self._interp('CFTBL3', X)
    def CFTBL4(self, X): return self._interp('CFTBL4', X)
    
    # TVAP / VPTBL aliases (Pressure -> Temp)
    def TVAP(self, X):   return self._interp('VPTBL', X) 
    def VPTBL(self, X):  return self._interp('VPTBL', X)
    
    # TBLVP (Temp -> Pressure)
    def TBLVP(self, X):  return self._interp('TBLVP', X)
    
    def ZTBLA(self, X):  return self._interp('ZTBLA', X)
    def ZTBLAP(self, X): return self._interp('ZTBLAP', X)
    def ZTBLD(self, X):  return self._interp('ZTBLD', X)
    def SHTBL1(self, X): return self._interp('SHTBL1', X)
    def SHTBL2(self, X): return self._interp('SHTBL2', X)
    def SHTBL3(self, X): return self._interp('SHTBL3', X)
    def SHTBL4(self, X): return self._interp('SHTBL4', X)
    def DELHV(self, X):  return self._interp('DELHV', X)
    def DELHL(self, X):  return self._interp('DELHL', X)

 
if __name__ == '__main__':
    a = PropertyTables()