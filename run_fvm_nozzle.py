"""Drive the axisymmetric FVM solver for the thruster sized by thruster_sizing.py.

Edit the input block below and run, or override from the command line::

    python run_fvm_nozzle.py --ni 240 --nj 90 --iters 30000
    python run_fvm_nozzle.py --euler            # inviscid, for comparison
    python run_fvm_nozzle.py --wall-temp 800    # cooled wall instead of adiabatic
    python run_fvm_nozzle.py --restart out/nozzle.npz --iters 10000
"""
import argparse
import os

import numpy as np
import metronos as mt

from fvm import (PerfectGas, LOX_LH2_OF5, ConicalNozzle, BellNozzle, Grid,
                 BoundaryConditions, NozzleSolver, quasi1d, post)

# ── Inputs ───────────────────────────────────────────────────────────────────
d_c      = mt.Quantity(8.3,   'mm')    # Chamber diameter
d_t      = mt.Quantity(0.29,  'mm')    # Throat diameter
expan    = 100.0                       # Area expansion ratio (A_e/A_t)
L_c      = mt.Quantity(10.0,  'mm')    # Chamber length (constant-radius section)
P_c      = mt.Quantity(8.45,  'bar')   # Chamber total pressure (from thruster_sizing.py)
P_a      = mt.Quantity(5.0,   'Torr')  # Ambient pressure
T_0      = mt.Quantity(3250., 'K')     # Chamber total temperature

theta_conv = 35.0                      # Convergent cone half-angle [deg]
theta_div  = 15.0                      # Divergent cone half-angle [deg]

# Combustion gas properties (LOX/LH2, O/F = 5.0), frozen and calorically perfect
GAS_KW = dict(LOX_LH2_OF5)


def build(args):
    gas = PerfectGas(**GAS_KW)
    Lc = L_c.to('m').value if args.l_chamber is None else args.l_chamber * 1e-3

    if args.bell:
        contour = BellNozzle(d_c.to('m').value / 2, d_t.to('m').value / 2,
                             expan, Lc, theta_conv=theta_conv)
    else:
        contour = ConicalNozzle(d_c.to('m').value / 2, d_t.to('m').value / 2,
                                expan, Lc,
                                theta_conv=theta_conv, theta_div=theta_div)

    grid = Grid.from_contour(contour, ni=args.ni, nj=args.nj,
                             wall_spacing=args.wall_spacing)
    wall = 'adiabatic' if args.wall_temp is None else args.wall_temp
    bcs = BoundaryConditions(p0=P_c.to('Pa').value, T0=T_0.to('K').value,
                             p_amb=P_a.to('Pa').value, wall=wall,
                             inlet=args.inlet)
    return gas, contour, grid, bcs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ni', type=int, default=200, help='axial cells')
    ap.add_argument('--nj', type=int, default=80, help='radial cells')
    ap.add_argument('--wall-spacing', type=float, default=0.004,
                    help='first cell height off the wall, as a fraction of local radius')
    ap.add_argument('--iters', type=int, default=20000)
    ap.add_argument('--tol', type=float, default=1e-7)
    ap.add_argument('--cfl', type=float, default=1.5)
    ap.add_argument('--flux', default='roe', choices=['roe', 'hllc', 'ausm'])
    ap.add_argument('--limiter', default='vanalbada',
                    choices=['vanalbada', 'minmod', 'vanleer', 'mc', 'none'])
    ap.add_argument('--integrator', default='rk5', choices=['rk1', 'rk3', 'rk5'])
    ap.add_argument('--inlet', default='stagnation',
                    choices=['stagnation', 'pressure', 'riemann'])
    ap.add_argument('--euler', action='store_true', help='inviscid run')
    ap.add_argument('--bell', action='store_true', help='bell contour instead of conical')
    ap.add_argument('--l-chamber', type=float, default=None, metavar='MM',
                    help='override the modelled chamber length [mm]. The chamber '
                         'is only a stagnation reservoir here, and at CR=819 it '
                         'sits at ~1 m/s, so its convective time scale is ~800x '
                         'the nozzle-s. Shortening it speeds convergence of the '
                         'global mass balance with almost no effect on thrust.')
    ap.add_argument('--wall-temp', type=float, default=None,
                    help='isothermal wall temperature [K] (default: adiabatic)')
    ap.add_argument('--out', default='out/nozzle')
    ap.add_argument('--restart', default=None)
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--print-every', type=int, default=500)
    args = ap.parse_args()

    np.seterr(all='ignore')
    gas, contour, grid, bcs = build(args)

    print(contour)
    print(grid.summary())
    print(f"\n{gas}")
    print(f"p0 = {P_c.to('bar').value:.3f} bar    T0 = {T_0.to('K').value:.1f} K    "
          f"p_amb = {P_a.to('Torr').value:.2f} Torr")
    print(f"contraction ratio = {contour.contraction_ratio:.1f}    "
          f"expansion ratio = {contour.area_ratio:.1f}")

    ideal = quasi1d.ideal_performance(contour, gas, bcs.p0, bcs.T0, bcs.p_amb)
    print(f"\nQuasi-1-D ideal: M_e = {ideal['M_exit']:.4f},  "
          f"mdot = {ideal['mdot'] * 1e3:.5f} g/s,  F = {ideal['thrust'] * 1e3:.3f} mN,  "
          f"Isp = {ideal['Isp']:.2f} s\n")

    solver = NozzleSolver(grid, gas, bcs, flux=args.flux, limiter=args.limiter,
                          viscous=not args.euler, cfl=args.cfl,
                          integrator=args.integrator)
    if args.restart and os.path.exists(args.restart):
        solver.initialize(quasi1d.initial_field(grid, gas, bcs.p0, bcs.T0, bcs.p_amb))
        solver.load(args.restart)
        print(f"restarted from {args.restart} at iteration {solver.iter}")
    else:
        solver.initialize(quasi1d.initial_field(grid, gas, bcs.p0, bcs.T0, bcs.p_amb))

    solver.run(max_iter=args.iters, tol=args.tol, print_every=args.print_every)

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    solver.save(args.out + '.npz')
    print(f"\ncheckpoint written to {args.out}.npz")

    post.report(solver, ideal)

    f = solver.fields()
    print(f"\nExit plane: M on axis {f['M'][-1, 0]:.3f},  "
          f"M at first cell off wall {f['M'][-1, -1]:.3f}")
    print(f"Boundary layer at exit: {post.boundary_layer_thickness(solver)[-1] * 1e3:.4f} mm "
          f"of a {contour.r_exit * 1e3:.4f} mm exit radius")

    paths = []
    if not args.no_plots:
        from fvm import plots
        paths = plots.save_all(solver, prefix=args.out, ideal=ideal)
        print("\nfigures written:")
        for p in paths:
            print("  " + p)

    md_path = args.out + '_RESULTS.md'
    with open(md_path, 'w', encoding='utf-8') as fh:
        fh.write(post.markdown_report(
            solver, ideal, figures=[os.path.basename(p) for p in paths],
            title='Axisymmetric Navier-Stokes solution'))
    print(f"\nsummary written to {md_path}")


if __name__ == '__main__':
    main()
