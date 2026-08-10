"""Post-processing: performance integrals, wall quantities and plots."""
import os

import numpy as np

from . import bc
from .thermo import G0

TINY = 1e-30
TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# integral quantities
# ---------------------------------------------------------------------------
def _iface_states(solver):
    """Primitive state on every i-face, from the ghosted cell array."""
    ng, ni, nj = solver.ng, solver.grid.ni, solver.grid.nj
    W = solver._fill(solver.U)
    js = slice(ng, ng + nj)
    return 0.5 * (W[:, ng - 1:ng + ni, js] + W[:, ng:ng + ni + 1, js])


def station_integrals(solver):
    """Mass flow and axial momentum+pressure flux at every i-face.

    Returns a dict of arrays of length ``ni+1``. The mass-flow array doubling
    as a conservation check is the cheapest useful convergence diagnostic
    there is: a converged nozzle solution has it flat to within a fraction of
    a percent.
    """
    g = solver.grid
    rho, u, v, p = _iface_states(solver)
    qn = u * g.nx_i + v * g.nr_i
    dA = TWO_PI * g.S_i

    mdot = np.sum(rho * qn * dA, axis=1)
    thrust = np.sum((rho * qn * u + (p - solver.bcs.p_amb) * g.nx_i) * dA, axis=1)
    area = np.sum(dA, axis=1)
    p_avg = np.sum(p * dA, axis=1) / np.maximum(area, TINY)

    # Mass-flow-weighted Mach and temperature: the honest one-dimensional
    # equivalent of the 2-D field, and the right thing to compare against the
    # quasi-1-D solution.
    a = np.sqrt(solver.gas.gamma * p / rho)
    M = np.hypot(u, v) / a
    T = p / (rho * solver.gas.R)
    w = rho * qn * dA
    wsum = np.sum(w, axis=1)
    safe = np.where(np.abs(wsum) < TINY, TINY, wsum)
    return dict(x=g.xf_i[:, 0], mdot=mdot, thrust=thrust, area=area, p_avg=p_avg,
                M_avg=np.sum(M * w, axis=1) / safe,
                T_avg=np.sum(T * w, axis=1) / safe)


def wall_quantities(solver):
    """Wall pressure, temperature, shear stress, heat flux and y+."""
    g, gas, bcs = solver.grid, solver.gas, solver.bcs
    ng, ni, nj = solver.ng, g.ni, g.nj
    W = solver._fill(solver.U)

    rho_c = W[0, ng:ng + ni, ng + nj - 1]
    u_c = W[1, ng:ng + ni, ng + nj - 1]
    v_c = W[2, ng:ng + ni, ng + nj - 1]
    p_c = W[3, ng:ng + ni, ng + nj - 1]
    T_c = p_c / (rho_c * gas.R)

    nx, nr = g.nx_j[:, -1], g.nr_j[:, -1]
    d = np.abs((g.xf_j[:, -1] - g.xc[:, -1]) * nx
               + (g.rf_j[:, -1] - g.rc[:, -1]) * nr)
    d = np.maximum(d, TINY)

    T_w = T_c if bcs.wall_is_adiabatic else np.full_like(T_c, bcs.T_wall)
    p_w = p_c
    rho_w = p_w / (gas.R * T_w)
    mu_w = gas.viscosity(T_w)
    kap_w = gas.conductivity(T_w)

    # Velocity gradient at the wall is purely wall-normal (no slip everywhere)
    dudn = -u_c / d
    dvdn = -v_c / d
    gux, gur = dudn * nx, dudn * nr
    gvx, gvr = dvdn * nx, dvdn * nr
    divV = gux + gvr                      # v = 0 on the wall, so the v/r term drops
    lam = -(2.0 / 3.0) * mu_w
    txx = 2.0 * mu_w * gux + lam * divV
    trr = 2.0 * mu_w * gvr + lam * divV
    txr = mu_w * (gur + gvx)
    tx = txx * nx + txr * nr
    tr = txr * nx + trr * nr
    tn = tx * nx + tr * nr
    tau_w = np.hypot(tx - tn * nx, tr - tn * nr)

    dTdn = np.zeros_like(T_c) if bcs.wall_is_adiabatic else (bcs.T_wall - T_c) / d
    q_w = -kap_w * dTdn                   # positive = heat flowing into the wall

    u_tau = np.sqrt(np.maximum(tau_w, 0.0) / rho_w)
    y_plus = rho_w * u_tau * d / mu_w

    return dict(x=g.xf_j[:, -1], r=g.rf_j[:, -1], p=p_w, T=T_w, tau=tau_w,
                q=q_w, y_plus=y_plus, u_tau=u_tau, d_wall=d,
                nx=nx, nr=nr, tx=tx, tr=tr, dA=TWO_PI * g.S_j[:, -1])


def performance(solver, ideal=None):
    """Thrust, mass flow, Isp, c*, Cf and the efficiencies against ideal."""
    g, gas, bcs = solver.grid, solver.gas, solver.bcs
    st = station_integrals(solver)
    wall = wall_quantities(solver)

    At = np.pi * g.contour.r_throat ** 2
    Ae = np.pi * g.contour.r_exit ** 2

    F = float(st["thrust"][-1])
    mdot = float(st["mdot"][-1])
    mdot_in = float(st["mdot"][0])

    # Mass conservation. The inlet plane is the slowest part of the domain to
    # settle (the chamber sits at M ~ 1e-3), so the throat-to-exit spread is
    # the meaningful number for the nozzle itself; the inlet figure is
    # reported separately rather than folded in.
    it = int(np.argmin(np.abs(st["x"] - g.contour.x_throat)))
    md_noz = st["mdot"][it:]
    mdot_spread = float(md_noz.max() - md_noz.min()) / max(abs(mdot), TINY)

    # Axial force the fluid exerts on the wall: pressure acts along the
    # outward normal, viscous traction is tau . n. Referenced to ambient, so
    # that the constant p_amb integrates to zero around the closed surface.
    F_wall_p = float(np.sum((wall["p"] - bcs.p_amb) * wall["nx"] * wall["dA"]))
    F_wall_v = float(np.sum(wall["tx"] * wall["dA"]))

    # Steady axial momentum balance over the control volume bounded by the
    # inlet plane, the exit plane and the wall:
    #     I_exit - I_inlet + (F_wall_p - F_wall_v) = 0
    # This is an independent check on the discretisation: it uses the wall
    # tractions, which never enter the exit-plane thrust integral.
    I_exit = F
    I_inlet = float(st["thrust"][0])
    imbalance = I_exit - I_inlet + (F_wall_p - F_wall_v)
    mom_residual = abs(imbalance) / max(abs(I_exit), TINY)

    out = dict(
        thrust=F, mdot=mdot, mdot_inlet=mdot_in,
        mdot_error=abs(mdot - mdot_in) / max(abs(mdot), TINY),
        mdot_spread_nozzle=mdot_spread,
        Isp=F / max(mdot * G0, TINY),
        cstar=bcs.p0 * At / max(mdot, TINY),
        Cf=F / max(bcs.p0 * At, TINY),
        throat_area=At, exit_area=Ae, area_ratio=Ae / At,
        wall_shear_axial=F_wall_v, wall_pressure_axial=F_wall_p,
        momentum_residual=mom_residual, momentum_imbalance=imbalance,
        heat_load=float(np.sum(wall["q"] * wall["dA"])),
        y_plus_max=float(np.nanmax(wall["y_plus"])),
    )

    if ideal is not None:
        out["thrust_ideal"] = ideal["thrust"]
        out["mdot_ideal"] = ideal["mdot"]
        out["Isp_ideal"] = ideal["Isp"]
        out["cstar_ideal"] = ideal["cstar"]
        out["Cd"] = mdot / max(ideal["mdot"], TINY)
        out["eta_thrust"] = F / max(ideal["thrust"], TINY)
        out["eta_Isp"] = out["Isp"] / max(ideal["Isp"], TINY)
        out["eta_cstar"] = out["cstar"] / max(ideal["cstar"], TINY)
    return out


def throat_reynolds(solver):
    """Reynolds number based on throat diameter and axis conditions there."""
    g, gas = solver.grid, solver.gas
    f = solver.fields()
    it = int(np.argmin(np.abs(g.xc[:, 0] - g.contour.x_throat)))
    rho = f["rho"][it, 0]
    u = f["speed"][it, 0]
    T = f["T"][it, 0]
    D = 2.0 * g.contour.r_throat
    return float(rho * u * D / gas.viscosity(T))


def boundary_layer_thickness(solver, frac=0.99):
    """Radial distance from the wall to ``frac`` of the peak speed, per station."""
    g = solver.grid
    f = solver.fields()
    q = f["speed"]
    qe = q.max(axis=1)
    delta = np.zeros(g.ni)
    for i in range(g.ni):
        prof = q[i]
        thr = frac * qe[i]
        idx = np.nonzero(prof >= thr)[0]
        if idx.size == 0:
            delta[i] = g.contour.r_wall(g.xc[i, 0])
        else:
            j = idx[-1]
            delta[i] = g.r_n[i, -1] - g.rc[i, j]
    return delta


def report(solver, ideal=None, quiet=False):
    """Formatted performance summary, in the style of thruster_sizing.py."""
    p = performance(solver, ideal)
    Re = throat_reynolds(solver)
    g = solver.grid

    L = []
    L.append("")
    L.append("FVM Nozzle Solution -- Performance Summary")
    L.append("-" * 78)
    L.append(f"{'Quantity':<28}{'CFD':>14}{'Ideal 1-D':>14}{'Ratio':>12}  Unit")
    L.append("-" * 78)

    def row(name, cfd, ide, unit, sc=1.0, fmt="{:>14.5g}"):
        a = fmt.format(cfd * sc)
        if ide is None:
            L.append(f"{name:<28}{a}{'-':>14}{'-':>12}  {unit}")
        else:
            b = fmt.format(ide * sc)
            L.append(f"{name:<28}{a}{b}{cfd / ide if ide else float('nan'):>12.4f}  {unit}")

    ig = ideal or {}
    row("Thrust", p["thrust"], ig.get("thrust"), "mN", 1e3)
    row("Mass flow", p["mdot"], ig.get("mdot"), "g/s", 1e3)
    row("Specific impulse", p["Isp"], ig.get("Isp"), "s")
    row("Characteristic velocity", p["cstar"], ig.get("cstar"), "m/s")
    row("Thrust coefficient", p["Cf"], ig.get("Cf"), "-")
    L.append("-" * 78)
    L.append(f"{'Throat Reynolds number':<28}{Re:>14.5g}")
    L.append(f"{'Max wall y+':<28}{p['y_plus_max']:>14.4g}")
    L.append(f"{'Wall axial pressure force':<28}{p['wall_pressure_axial'] * 1e3:>14.5g}  mN")
    L.append(f"{'Wall axial viscous force':<28}{p['wall_shear_axial'] * 1e3:>14.5g}  mN")
    L.append(f"{'Integrated wall heat load':<28}{p['heat_load']:>14.5g}  W")
    L.append(f"{'Mass spread, throat->exit':<28}{p['mdot_spread_nozzle'] * 100:>14.4g}  %")
    L.append(f"{'Mass error, inlet vs exit':<28}{p['mdot_error'] * 100:>14.4g}  %")
    L.append(f"{'Momentum balance residual':<28}{p['momentum_residual'] * 100:>14.4g}  %")
    L.append("-" * 78)
    if solver.history:
        L.append(f"{'Iterations':<28}{solver.iter:>14d}")
        L.append(f"{'Residual drop':<28}{solver.history[-1][2]:>14.3e}")
    L.append(f"{'Grid':<28}{f'{g.ni} x {g.nj}':>14}")
    L.append("-" * 78)

    text = "\n".join(L)
    if not quiet:
        print(text)
    return text


def markdown_report(solver, ideal=None, figures=(), title="FVM nozzle results"):
    """A shareable summary of the run, with the figures embedded."""
    p = performance(solver, ideal)
    g = solver.grid
    b = solver.bcs
    Re = throat_reynolds(solver)
    delta = boundary_layer_thickness(solver)[-1]
    f = solver.fields()

    L = [f"# {title}", ""]
    L.append("## Case")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    L.append(f"| Throat diameter | {2e3 * g.contour.r_throat:.4f} mm |")
    L.append(f"| Expansion ratio | {g.contour.area_ratio:.1f} |")
    L.append(f"| Contraction ratio | {g.contour.contraction_ratio:.1f} |")
    L.append(f"| Chamber total pressure | {b.p0 / 1e5:.3f} bar |")
    L.append(f"| Chamber total temperature | {b.T0:.0f} K |")
    L.append(f"| Ambient pressure | {b.p_amb / 133.322:.2f} Torr |")
    L.append(f"| Wall | {'adiabatic' if b.wall_is_adiabatic else f'isothermal, {b.T_wall:.0f} K'} |")
    L.append(f"| Grid | {g.ni} x {g.nj} = {g.ni * g.nj} cells |")
    L.append(f"| Flux / limiter | {solver.flux_name} / {solver.limiter_name} |")
    L.append(f"| Physics | {'Navier-Stokes (laminar)' if solver.viscous else 'Euler'} |")
    L.append("")

    L.append("## Performance")
    L.append("")
    L.append("| Quantity | CFD | Ideal 1-D | Ratio |")
    L.append("|---|---|---|---|")

    def row(name, cfd, key, unit, sc=1.0, fmt="{:.5g}"):
        ide = (ideal or {}).get(key)
        a = fmt.format(cfd * sc) + f" {unit}"
        if ide is None:
            L.append(f"| {name} | {a} | - | - |")
        else:
            L.append(f"| {name} | {a} | {fmt.format(ide * sc)} {unit} | "
                     f"{cfd / ide:.4f} |")

    row("Thrust", p["thrust"], "thrust", "mN", 1e3)
    row("Mass flow", p["mdot"], "mdot", "g/s", 1e3)
    row("Specific impulse", p["Isp"], "Isp", "s")
    row("Characteristic velocity", p["cstar"], "cstar", "m/s")
    row("Thrust coefficient", p["Cf"], "Cf", "-")
    L.append("")

    L.append("## Viscous diagnostics")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    L.append(f"| Throat Reynolds number | {Re:.0f} |")
    L.append(f"| Exit Mach, axis | {f['M'][-1, 0]:.3f} |")
    L.append(f"| Exit Mach, first cell off wall | {f['M'][-1, -1]:.3f} |")
    L.append(f"| Boundary layer at exit | {delta * 1e3:.4f} mm of "
             f"{g.contour.r_exit * 1e3:.4f} mm exit radius |")
    L.append(f"| Max wall y+ | {p['y_plus_max']:.3f} |")
    L.append(f"| Wall axial viscous force | {p['wall_shear_axial'] * 1e3:.4g} mN |")
    if not b.wall_is_adiabatic:
        L.append(f"| Integrated wall heat load | {p['heat_load']:.4g} W |")
    L.append("")

    L.append("## Numerical quality")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    L.append(f"| Iterations | {solver.iter} |")
    if solver.history:
        L.append(f"| Residual drop | {solver.history[-1][2]:.3e} |")
    L.append(f"| Mass spread, throat to exit | {p['mdot_spread_nozzle'] * 100:.4g} % |")
    L.append(f"| Mass error, inlet vs exit | {p['mdot_error'] * 100:.4g} % |")
    L.append(f"| Momentum balance residual | {p['momentum_residual'] * 100:.4g} % |")
    L.append("")

    if figures:
        L.append("## Figures")
        L.append("")
        for path in figures:
            name = os.path.basename(path).rsplit(".", 1)[0]
            L.append(f"### {name.replace('_', ' ')}")
            L.append("")
            L.append(f"![{name}]({path.replace(os.sep, '/')})")
            L.append("")
    return "\n".join(L)
