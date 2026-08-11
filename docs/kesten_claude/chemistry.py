import json
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from properties import PropertyTables

@dataclass
class Species:
    name: str
    mw: float
    diffusivity: float = 0.0 
    thermo_config: Dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            name=data['name'],
            mw=data['molecular_weight'],
            diffusivity=data.get('transport', {}).get('diffusivity', 0.0),
            thermo_config=data.get('thermo', {})
        )

    def get_cp(self, T: float, props: PropertyTables) -> float:
        """Looks up specific heat using the table name defined in JSON."""
        table_name = self.thermo_config.get('cp_table')
        if not table_name: return 0.0
        # Dynamically call the method on the props object (e.g., props.CFTBL1(T))
        method = getattr(props, table_name, None)
        if method:
            return float(method(T))
        return 0.0

    def get_enthalpy_reaction(self, T: float, props: PropertyTables) -> float:
        """Legacy support for heat of reaction tables (H3TBL/H4TBL)."""
        table_name = self.thermo_config.get('h_reaction_table')
        if not table_name: return 0.0
        method = getattr(props, table_name, None)
        if method:
            return float(method(T))
        return 0.0


@dataclass
class Reaction:
    id: str
    reactants: Dict[str, float]
    products: Dict[str, float]
    A: float
    b: float
    Ea_R: float
    orders: Dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict):
        kin = data.get('kinetics', {})
        return cls(
            id=data['id'],
            reactants=data.get('reactants', {}),
            products=data.get('products', {}),
            A=kin.get('A', 0.0),
            b=kin.get('b', 0.0),
            Ea_R=kin.get('Ea_R', 0.0),
            orders=data.get('orders', {})
        )

    def get_net_stoichiometry(self, species_name: str) -> float:
        return self.products.get(species_name, 0.0) - self.reactants.get(species_name, 0.0)

    def get_rate(self, T: float, concentrations: Dict[str, float]) -> float:
        if T <= 0: return 0.0
        k = self.A * (T ** self.b) * np.exp(-self.Ea_R / T)
        
        # Law of Mass Action (or explicit orders)
        drivers = self.orders if self.orders else self.reactants
        rate_prod = 1.0
        for sp, order in drivers.items():
            conc = max(concentrations.get(sp, 0.0), 1e-30)
            rate_prod *= (conc ** order)
        return k * rate_prod


class Mechanism:
    def __init__(self, json_path: str):
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        self.species = [Species.from_dict(s) for s in data['species']]
        self.reactions = [Reaction.from_dict(r) for r in data['reactions']]
        
        self.species_map = {s.name: s for s in self.species}
        self.reaction_map = {r.id: r for r in self.reactions}

    def get_species(self, name):
        return self.species_map.get(name)

    def get_reaction(self, id):
        return self.reaction_map.get(id)