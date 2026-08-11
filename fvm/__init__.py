"""Axisymmetric finite-volume Navier-Stokes solver for small nozzles.

A second-order cell-centred FVM on a structured, body-fitted grid, written to
resolve the flowfield inside the micro-thrusters that ``thruster_sizing.py``
sizes algebraically. At a 0.29 mm throat the boundary layer occupies a large
fraction of the passage, so the viscous terms are not a correction -- they set
the answer.

Typical use::

    from fvm import (PerfectGas, LOX_LH2_OF5, ConicalNozzle, Grid,
                     BoundaryConditions, NozzleSolver, quasi1d, post)

    gas = PerfectGas(**LOX_LH2_OF5)
    contour = ConicalNozzle(r_chamber=4.15e-3, r_throat=0.145e-3,
                            area_ratio=100.0, L_chamber=10e-3)
    grid = Grid.from_contour(contour, ni=240, nj=80)
    bcs = BoundaryConditions(p0=6.9e5, T0=3250.0, p_amb=666.6)

    sol = NozzleSolver(grid, gas, bcs)
    sol.initialize(quasi1d.initial_field(grid, gas, bcs.p0, bcs.T0, bcs.p_amb))
    sol.run(max_iter=20000, tol=1e-6)
    print(post.report(sol))
"""

from .thermo import PerfectGas, LOX_LH2_OF5, COLD_N2, G0, R_UNIVERSAL
from .geometry import NozzleContour, ConicalNozzle, BellNozzle
from .grid import Grid
from .bc import BoundaryConditions
from .catbed import (CatalystBed, PackedSpheres, PackedCylinders,
                     ReticulatedFoam, mesh_to_diameter)
from .mechanism import HydrazineShell405, get_mechanism
from .plugflow import PlugFlowReactor, PlugFlowSolution
from .thruster import ThrusterSystem, ThrusterSolution, vapor_region_inlet
from .solver import NozzleSolver
from . import quasi1d, post, riemann, reconstruct, viscous, chem, catbed, mechanism, plugflow, thruster

__all__ = [
    "PerfectGas", "LOX_LH2_OF5", "COLD_N2", "G0", "R_UNIVERSAL",
    "NozzleContour", "ConicalNozzle", "BellNozzle",
    "Grid", "BoundaryConditions", "NozzleSolver",
    "CatalystBed", "PackedSpheres", "PackedCylinders",
    "ReticulatedFoam", "mesh_to_diameter",
    "HydrazineShell405", "get_mechanism",
    "PlugFlowReactor", "PlugFlowSolution",
    "ThrusterSystem", "ThrusterSolution", "vapor_region_inlet",
    "quasi1d", "post", "riemann", "reconstruct", "viscous",
    "chem", "catbed", "mechanism", "plugflow", "thruster",
]

__version__ = "0.1.0"
