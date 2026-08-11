from dataclasses import dataclass
from properties import PropertyTables
# Ensure chemistry.py is in the same folder
from chemistry import Mechanism 

@dataclass
class ReactorConfig:
    # --- Integration Controls ---
    compare: bool = True
    print_rows: bool = True
    refine_event: bool = True
    max_steps: int = 8000
    
    # --- Geometry ---
    G0: float = 1.51
    Z0: float = 0.0
    ZEND: float = 0.25
    FC: float = 0.0
    GATZ0: float = 3.00
    
    # --- BCs ---
    HF: float = 0.0
    TF: float = 530.0
    PRES: float = 111.4
    
    # --- Constants ---
    R: float = 10.73
    CFL: float = 0.7332
    KP: float = 0.4e-4
    C1: float = 1.0
    
    # --- Mechanism & Chemistry ---
    # This is the new field main.py is looking for
    mechanism_file: str = "kinetics.json"
    mechanism: Mechanism = None
    
    # --- Legacy Bridges (Populated from Mechanism) ---
    # These allow the legacy LiquidPhase logic (SLOPE/PARAM) to keep working
    WM1: float = 0.0; WM2: float = 0.0; WM3: float = 0.0; WM4: float = 0.0
    DIF3: float = 0.0; DIF4: float = 0.0
    ALPHA1: float = 0.0; ALPHA2: float = 0.0
    ALPHA3: float = 2.14e10
    AGM: float = 0.0; BGM: float = 0.0
    CGM: float = 33000.0
    EN1: float = 1.0; EN3: float = -1.6
    ENMX1: float = 200.0; ENMX2: float = 40.0; ENMX3: float = 80.0
    
    # --- Calculated Properties ---
    TVAP: float = 0.0
    HL: float = 0.0
    HV: float = 0.0

    @classmethod
    def from_defaults(cls, props: PropertyTables):
        cfg = cls()
        
        # 1. Load Mechanism
        # This will fail if chemistry.py or kinetics.json are missing
        cfg.mechanism = Mechanism(cfg.mechanism_file)
        
        # 2. Map Chemistry to Legacy Fields (Bridging for SLOPE/PARAM)
        # We map specific species to the legacy "1, 2, 3, 4" slots
        hydrazine = cfg.mechanism.get_species("N2H4")
        ammonia = cfg.mechanism.get_species("NH3")
        nitrogen = cfg.mechanism.get_species("N2")
        hydrogen = cfg.mechanism.get_species("H2")
        
        if hydrazine: cfg.WM4, cfg.DIF4 = hydrazine.mw, hydrazine.diffusivity
        if ammonia:   cfg.WM3, cfg.DIF3 = ammonia.mw, ammonia.diffusivity
        if nitrogen:  cfg.WM2 = nitrogen.mw
        if hydrogen:  cfg.WM1 = hydrogen.mw
        
        r1 = cfg.mechanism.get_reaction("R1")
        if r1: cfg.ALPHA1, cfg.AGM = r1.A, r1.Ea_R
        
        r2 = cfg.mechanism.get_reaction("R2")
        if r2: 
            cfg.ALPHA2, cfg.BGM = r2.A, r2.Ea_R
            cfg.EN3 = r2.orders.get("H2", -1.6)

        # 3. Standard Calcs
        cfg.TVAP = props.TVAP(cfg.PRES)
        cfg.HL = (cfg.TVAP - cfg.TF) * cfg.CFL
        cfg.HV = cfg.HL + props.DHVST(cfg.TVAP) - props.DHLVST(cfg.TVAP)
        cfg.GATZ0 = cfg.G0 + cfg.FC * cfg.Z0
        
        return cfg
