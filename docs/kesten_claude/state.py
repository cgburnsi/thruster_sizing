import numpy as np
from dataclasses import dataclass, field
from typing import Optional

@dataclass
class ReactorState:
    """
    Holds the time-dependent (or spatially-dependent) state arrays 
    for the reactor simulation.
    """
    max_steps: int
    
    # State Arrays (Initialized to zeros or NaN)
    Z: np.ndarray = field(init=False)
    H: np.ndarray = field(init=False)
    TEMP: np.ndarray = field(init=False)
    DHDZ: np.ndarray = field(init=False)
    DERIV: np.ndarray = field(init=False)
    WFV: np.ndarray = field(init=False)
    PRES: np.ndarray = field(init=False)
    
    # Species Concentrations (Vapor Region)
    C1: np.ndarray = field(init=False) # H2
    C2: np.ndarray = field(init=False) # N2
    C3: np.ndarray = field(init=False) # NH3
    C4: np.ndarray = field(init=False) # N2H4
    
    # Gradients
    GRAD3: np.ndarray = field(init=False)
    TGRAD: np.ndarray = field(init=False)

    # Tracking Indices
    ii: int = 0
    ii_lqvp_end: Optional[int] = None
    ii_vap_end: Optional[int] = None
    
    # Event Markers
    z_lv: float = 0.0      # Z where liquid-vapor transition starts
    dz_last: float = 0.0   # Last dz step size (for handing off between regions)

    def __post_init__(self):
        """Allocate memory for arrays based on max_steps."""
        n = self.max_steps
        self.Z = np.zeros(n)
        self.H = np.zeros(n)
        self.TEMP = np.zeros(n)
        self.DHDZ = np.zeros(n)
        self.DERIV = np.zeros(n)
        self.WFV = np.zeros(n)
        self.PRES = np.full(n, np.nan)
        
        # Initialize species/gradients with NaN to avoid confusion in plots
        self.C1 = np.full(n, np.nan)
        self.C2 = np.full(n, np.nan)
        self.C3 = np.full(n, np.nan)
        self.C4 = np.full(n, np.nan)
        self.GRAD3 = np.full(n, np.nan)
        self.TGRAD = np.full(n, np.nan)
    
    @property
    def current_depth(self):
        return self.Z[self.ii]

    @property
    def current_enthalpy(self):
        return self.H[self.ii]

    def pack_diag(self, idx, diag_dict):
        """Helper to store diagnostic values into the arrays at index idx."""
        if "temp" in diag_dict: self.TEMP[idx] = diag_dict["temp"]
        if "deriv" in diag_dict: self.DERIV[idx] = diag_dict["deriv"]
        if "grad3" in diag_dict: self.GRAD3[idx] = diag_dict["grad3"]
        if "tgrad" in diag_dict: self.TGRAD[idx] = diag_dict["tgrad"]
