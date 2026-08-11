import json
import numpy as np
from pathlib import Path

class Species:
    def __init__(self, rec):
        self.R          = 8.314462618                   # [J/mol-K] Universal Gas Constant
        self.key        = rec['key']                    # [-] Identification of species
        self.name       = rec.get('name', self.key)     # [-] Official name for the species
        self.aliases    = rec.get('aliases', [])        # [-] Other names for the species
        self.formula    = rec.get('formula')            # [-] Chemical Formula for Species
        self.phase      = rec.get('phase', 'g')         # [-] phase (liquid, solid, gas)
        self.mw         = rec.get('mw')                 # [g/mol] Molecular Mass of Species
        self.hf298      = rec.get('hf298')              # [J/mol] Assigned Enthalpy at 298.15 K for stadard state
        self.hf298_h0   = rec.get('hf298_h0')           # [J/mol] 
        self.elements   = rec.get('elements', {})       # [-] Individual elements in the species
        
        # Process the record
        self._get_thermo_blocks(rec)        # Get the thermodynamics coefficients
        self._get_visc_blocks(rec)          # Get the viscosity blocks if present in the record
        self._get_cond_blocks(rec)          # Get the thermal conductivity blocks if present in the record

    # ----- Useful Helper Methods -------------------------------------------------------------------------------------
    @property
    def has_viscosity(self):        return bool(self.visc_blocks)
    
    @property
    def has_conductivity(self):     return bool(self.cond_blocks)
    
    @property
    def T_min(self):                return self.blocks[0]['T_low']
    
    @property
    def T_max(self):                return self.blocks[-1]['T_high']
    
    # ----- Species Class PRIVATE METHODS ----------------------------------------------------------------------------- 
    def _pick_block(self, blocks, T, label):
        if not blocks:  raise ValueError('%s: No Coefficient Blocks Defined' % self.key)
        if T <= 0.0:    raise ValueError('%s: T must be > 0 K (got %g)' % (label, T)) 
        
        minT = self.blocks[0]['T_low']
        maxT = self.blocks[-1]['T_high']
        if T < minT or T > maxT: raise ValueError('%s: T=%g K is outside valid range [%g, %g] K' % (label, T, minT, maxT))
        
        for i, block in enumerate(self.blocks):
            low, high = block['T_low'], block['T_high']
            # include the high endpoint only for the final block
            if (low <=T < high) or (i == len(self.blocks) -1 and T == high):
                return block
        
        # if you get here the blocks in the data are incorrect
        raise ValueError('%s: No NASA-9 block matched T=%g K; check your ranges.' % (self.key, T))

    def _coeffs_for_T(self, T):
        return self._pick_block(self.blocks, T, self.key)
    
    def _get_thermo_blocks(self, rec):
        ''' Get the thermo property coefficients from the record that was loaded. '''
        blocks = rec.get('thermo', [])
        if not blocks: raise ValueError('%s: missing "thermo" blocks' % self.key)
        
        for block in blocks:
            if len(block.get('a', [])) != 7: raise ValueError('%s: NASA-9 block must have 7 "a" coefficients' % self.key)
            if len(block.get('b', [])) != 2: raise ValueError('%s: NASA-9 block must have 2 "b" coefficients' % self.key)
        self.blocks = blocks
        
    def _get_visc_blocks(self, rec):
        ''' Get dynamics viscosity coefficients from the record that was loaded, if available. '''
        blocks = rec.get('visc', [])
        if not blocks: 
            self.visc_blocks = None
            return
        self.visc_blocks = rec.get('visc') # Get the coefficient blocks if you get to this line.
    
    def _get_cond_blocks(self, rec):
        self.cond_blocks = rec.get('cond') or None   # None if missing/empty

    # ----- Species Class PUBLIC API ----------------------------------------------------------------------------------
    def Cp(self, T):
        block                       = self._coeffs_for_T(T)
        a1, a2, a3, a4, a5, a6, a7  = block['a']
        b1, b2                      = block['b']
        
        term    = np.zeros(7)
        term[0] = a1 / T**2
        term[1] = a2 / T
        term[2] = a3
        term[3] = a4 * T
        term[4] = a5 * T**2
        term[5] = a6 * T**3
        term[6] = a7 * T**4

        return np.sum(term) * self.R

    def H(self, T):
        block                       = self._coeffs_for_T(T)
        a1, a2, a3, a4, a5, a6, a7  = block['a']
        b1, b2                      = block['b']
        
        term    = np.zeros(8)
        term[0] = -a1/T**2
        term[1] =  a2 * np.log(T) / T
        term[2] =  a3
        term[3] =  a4 * T / 2.0
        term[4] =  a5 * T**2 / 3.0
        term[5] =  a6 * T**3 / 4.0
        term[6] =  a7 * T**4 / 5.0
        term[7] =  b1 / T
        
        return np.sum(term) * self.R * T

    def S(self, T):
        block                       = self._coeffs_for_T(T)
        a1, a2, a3, a4, a5, a6, a7  = block['a']
        b1, b2                      = block['b']
        
        term    = np.zeros(8)
        term[0] = -a1 / (2.0 * T**2)
        term[1] = -a2 / T
        term[2] =  a3 * np.log(T)
        term[3] =  a4 * T 
        term[4] =  a5 * T**2 / 2.0
        term[5] =  a6 * T**3 / 3.0
        term[6] =  a7 * T**4 / 4.0
        term[7] =  b2
        
        return np.sum(term) * self.R
        
    def mu(self, T, default=None):
        """Dynamic viscosity [Pa·s]; returns `default` if no data or block malformed."""
        blocks = getattr(self, 'visc_blocks', None)
        if not blocks:
            return default
        try:
            blk = self._pick_block(blocks, T, self.key + " (visc)")
        except ValueError:
            return default
    
        coeffs = blk.get('mu') or blk.get('coeffs') or blk.get('a')
        if not coeffs or len(coeffs) != 4:
            # optional: print(f"{self.key} visc block missing 4 coeffs; keys={list(blk.keys())}")
            return default
    
        A, B, C, D = coeffs
        form = blk.get('form', 'mu_log10_microP')
        if form == 'mu_log10_microP':
            # μ[Pa·s] = 1e-7 * 10^(A + B/T + C/T^2 + D)
            return 1e-7 * (10.0 ** (A + B/T + C/(T*T) + D))
        elif form == 'mu_log10_microP_logT':
            return 1e-7 * (10.0 ** (A + B/T + C/(T*T) + D*np.log10(T)))
        # unknown form → be forgiving
        return default
    
    def k(self, T, default=None):
        """Thermal conductivity [W/(m·K)]; returns `default` if no data or block malformed."""
        blocks = getattr(self, 'cond_blocks', None)
        if not blocks:
            return default
        try:
            blk = self._pick_block(blocks, T, self.key + " (cond)")
        except ValueError:
            return default
    
        coeffs = blk.get('k') or blk.get('coeffs') or blk.get('a')
        if not coeffs or len(coeffs) != 4:
            return default
    
        A, B, C, D = coeffs
        form = blk.get('form', 'k_ln_WmK')
        if form == 'k_ln_WmK':
            return float(np.exp(A + B/T + C/(T*T) + D))
        elif form == 'k_ln_mWmK':
            return 1e-3 * float(np.exp(A + B/T + C/(T*T) + D))
        return default
   

        



class SpeciesLibrary:
    ''' In-memory collection of multiple species with key/alias lookup and simple loaders. '''
    def __init__(self):
        self.by_key   = {}
        self.by_alias = {}
    
    def __iter__(self):         return iter(self.by_key.values())    
    def __len__(self):          return len(self.by_key)    
    def __contains__(self, k):  
        k = str(k)
        return (k in self.by_key) or (k.lower() in self.by_alias)
    
    # ----- SpeciesLibrary Class PUBLIC API ---------------------------------------------------------------------------
    def add(self, sp_or_rec, overwrite = False):
        """Add a Species or a raw record. First-alias-wins policy (setdefault)."""
        sp = sp_or_rec if isinstance(sp_or_rec, Species) else Species(sp_or_rec)

        if (not overwrite) and (sp.key in self.by_key):
            raise ValueError("Duplicate species key: %s" % sp.key)

        self.by_key[sp.key] = sp

        alias_list = list(sp.aliases) + [sp.key]
        if sp.formula:
            alias_list.append(sp.formula)
        for a in alias_list:
            if a:
                self.by_alias.setdefault(a.lower(), sp.key)  # keep first claimant

        return sp

    def get(self, key_or_alias, default=None):
        """Return Species by canonical key or alias (case-insensitive)."""
        sp = self.by_key.get(key_or_alias)
        if sp is not None:
            return sp
        k = self.by_alias.get(str(key_or_alias).lower())
        return self.by_key.get(k, default)

    # ---- convenience loaders (optional I/O) ----
    @staticmethod
    def _default_loader(path):
        p = Path(path)
        txt = p.read_text(encoding="utf-8-sig").strip()
        if not txt:
            raise ValueError("Empty file: %s" % p)
        return json.loads(txt)

    @classmethod
    def from_files(cls, *paths, loader=None):
        lib = cls()
        load = loader or cls._default_loader
        for p in paths:
            lib.add(load(p))
        return lib

    @classmethod
    def from_dir(cls, directory, pattern="*.json", loader=None):
        lib = cls()
        load = loader or cls._default_loader
        for fp in sorted(Path(directory).glob(pattern)):
            lib.add(load(fp))
        return lib
    
    

if __name__ == '__main__':
    
    # 1) Build the library from a directory of JSON files
    lib = SpeciesLibrary.from_dir(Path("species"))  # or Path(r"C:\full\path\species")
    print(f"Loaded {len(lib)} species:")
    print([sp.key for sp in lib])  # thanks to __iter__, this yields Species objects

    # 2) Lookups: canonical key or alias (case-insensitive)
    nitrogen = lib.get("N2")         # by key
    hydrogen = lib.get("hydrogen")       # by alias
    print("nitrogen key:", nitrogen.key, "| hydrogen key:", hydrogen.key)

    # 3) Iterate and compute Cp at a temperature (skip out-of-range cleanly)
    T = 298.15
    print(f"\nCp at {T} K:")
    for sp in lib:
        try:
            Cpi = sp.Cp(T)
            Hi  = sp.H(T)
            Si  = sp.S(T)
        except ValueError as e:  # raised by your strict _coeffs_for_T
            print(f"  skip {sp.key}: {e}")
            continue
        print(f"  {sp.key}: Cp = {Cpi:.1f} J/(mol·K)")
        print(f"  {sp.key}: H  = {Hi:.1f} J/mol")
        print(f"  {sp.key}: S  = {Si:.1f} J/mol")

    # 4) Membership checks include aliases (from our __contains__ tweak)
    print("\nMembership tests:")
    print("  'H2' in lib?        ", "H2" in lib)
    print("  'N2' in lib?", "N2" in lib)

    # 5) Safe get with default (mirrors dict.get)
    print("\nUnknown species returns default:")
    print("  lib.get('Nope', None) ->", lib.get("Nope", None))
    
    
    nh3 = lib.get("NH3")
    T = 600.0
    
    mu = nh3.mu(T)              # None if absent or out-of-range
    k  = nh3.k(T, default=np.nan)
    
    #print("μ =", mu if mu is not None else "N/A")
    #print("k =", k)



