# Thruster Sizing

A Python tool for sizing small chemical thrusters. Given a target thrust level and nozzle geometry, it solves for exit Mach number and chamber pressure, then computes key performance parameters.

## Features

- Newton-Raphson solver for supersonic exit Mach number from area expansion ratio
- Newton-Raphson solver for chamber pressure from thrust, geometry, and ambient conditions
- Outputs a formatted performance and design summary (SI and English units)
- Catalyst bed surface area and loading estimates
- Chamber wall thickness estimate from yield strength and max pressure

## Requirements

Install Python dependencies:

```
pip install -r requirements.txt
```

> **Note:** `metronos` is installed directly from GitHub (see `requirements.txt`). It is not available on PyPI.

## Usage

Edit the inputs at the top of the `__main__` block in `thruster_sizing.py`, then run:

```
python thruster_sizing.py
```

### Key Inputs

| Parameter | Description |
|---|---|
| `d_t` | Throat diameter |
| `expan` | Area expansion ratio (A_e / A_t) |
| `F_set` | Design thrust |
| `P_a` | Ambient pressure |
| `k` | Combustion gas specific heat ratio |
| `cstar_ideal` | Ideal characteristic velocity |
| `eta` | C* efficiency |

### Example Output

```
Thruster Performance and Design Summary
--------------------------------------------------------------------------------
Parameter                        SI Value Unit              Eng Value Unit
--------------------------------------------------------------------------------
Chamber Diameter, d_c              8.3000 mm                   0.3268 in
Throat Diameter, d_t               0.2900 mm                   0.0114 in
Exit Diameter, d_e                 2.9000 mm                   0.1142 in
...
```

## Propellant

The default combustion properties are for **LOX/LH2** at a mixture ratio of O/F = 5.0. These are derived from thermochemical equilibrium analysis (e.g., NASA CEA). To use a different propellant, update `k`, `cstar_ideal`, `MW`, and `T0` in the `__main__` block.

| Property | Value | Description |
|---|---|---|
| `k` | 1.26 | Specific heat ratio |
| `cstar_ideal` | 2350 m/s | Ideal characteristic velocity |
| `MW` | 11.8 g/mol | Combustion gas molecular weight |
| `T0` | 3250 K | Adiabatic flame temperature |
| `eta` | 0.95 | C* efficiency |
