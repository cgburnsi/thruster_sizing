"""Validation against Hayn's micronozzle experiments (NASA TM-77730).

Dieter Hayn, "Beitraege zur Leistungsermittlung von Mikroduesen", TU Munich
dissertation 1983, translated as NASA TM-77730 (1984). Twelve cold-gas (N2)
micronozzles were formed galvanoplastically over a turned core, then machined
down through five area-ratio steps and fired at five chamber pressures --
1200 individual tests.

Nozzle 2 at its epsilon = 100 step is, to within about 2%, the thruster that
``thruster_sizing.py`` sizes:

    throat radius   0.142 mm   (vs 0.145 mm)
    exit diameter   2.84 mm    (vs 2.90 mm)
    area ratio      99.65      (vs 100)
    divergence      15 deg     (vs 15 deg)

Geometry is from Table 4, p.107. The table gives throat radius R_T, exit
diameter D and an overall length l measured from a fixed datum. Differencing
l and D across the five epsilon steps recovers a 15 deg divergent half-angle
and puts the throat 8.45 mm downstream of the datum; with the 5.0 mm chamber
(Bk_L) that leaves 3.45 mm of convergent section, i.e. a ~28 deg convergent
half-angle from the 2.0 mm chamber radius (Bk_R).

Two caveats on what this validates:

* The Reynolds numbers are not the thruster's. Cold N2 at 293 K runs
  Re_throat = 8700 (2 bar) to 87000 (20 bar); the LOX/LH2 thruster runs 1160.
  Identical geometry, different viscous regime, so agreement here does not
  transfer directly to the hot case -- it validates the solver, not the
  thruster prediction.
* A laminar solver is nonetheless the right tool. Hayn addresses this
  explicitly (p.31): the strong favourable pressure gradient puts the Boldman
  acceleration parameter several times above critical, relaminarising the
  boundary layer, so all of his own boundary-layer calculations assume a
  laminar profile too.

Usage::

    python validate_hayn.py --pc 10 --iters 25000
    python validate_hayn.py --pc 2 --p-amb 1.0
"""
import argparse
import os

import numpy as np

from fvm import (PerfectGas, COLD_N2, ConicalNozzle, Grid,
                 BoundaryConditions, NozzleSolver, quasi1d, post)

# ── Hayn nozzle 2, epsilon = 100 step (Table 4, p.107) ───────────────────────
R_THROAT = 0.142e-3      # m   throat radius
AREA_RATIO = 99.65       # -   measured area ratio at this step
R_CHAMBER = 2.0e-3       # m   Bk_R
L_CHAMBER = 5.0e-3       # m   Bk_L
THETA_DIV = 15.0         # deg alpha_e
THETA_CONV = 28.0        # deg inferred from the length budget (see docstring)
T0 = 293.0               # K   cold gas, room temperature

#: Chamber pressures Hayn tested, in bar (p.110)
PC_TESTED = (2.0, 5.0, 10.0, 15.0, 20.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pc', type=float, default=10.0,
                    help=f'chamber pressure [bar]; Hayn tested {PC_TESTED}')
    ap.add_argument('--p-amb', type=float, default=1.0,
                    help='vacuum chamber pressure [Torr]; his facility ran 1e-2 to 10 Torr')
    ap.add_argument('--ni', type=int, default=200)
    ap.add_argument('--nj', type=int, default=96)
    ap.add_argument('--wall-spacing', type=float, default=0.0018)
    ap.add_argument('--iters', type=int, default=25000)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--cfl', type=float, default=1.5)
    ap.add_argument('--euler', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--print-every', type=int, default=2000)
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--checkpoint-every', type=int, default=1000,
                    help='write a restart file every N iterations (0 disables)')
    ap.add_argument('--fresh', action='store_true',
                    help='ignore any existing checkpoint and start over')
    args = ap.parse_args()

    np.seterr(all='ignore')
    out = args.out or f'out/hayn_pc{args.pc:g}bar'

    gas = PerfectGas(**COLD_N2)
    contour = ConicalNozzle(R_CHAMBER, R_THROAT, AREA_RATIO, L_CHAMBER,
                            theta_conv=THETA_CONV, theta_div=THETA_DIV)
    grid = Grid.from_contour(contour, ni=args.ni, nj=args.nj,
                             wall_spacing=args.wall_spacing)
    p0 = args.pc * 1e5
    p_amb = args.p_amb * 133.322
    bcs = BoundaryConditions(p0=p0, T0=T0, p_amb=p_amb, wall='adiabatic')

    print('Hayn nozzle 2, epsilon = 100 step   (NASA TM-77730, Table 4)')
    print(contour)
    print(grid.summary())
    print(f'\n{gas}')
    print(f'p0 = {args.pc:g} bar   T0 = {T0:.0f} K   '
          f'p_amb = {args.p_amb:g} Torr ({p_amb:.1f} Pa)')

    ideal = quasi1d.ideal_performance(contour, gas, p0, T0, p_amb)
    print(f'\nQuasi-1-D ideal: M_e = {ideal["M_exit"]:.3f}, '
          f'mdot = {ideal["mdot"] * 1e3:.5f} g/s, F = {ideal["thrust"] * 1e3:.3f} mN, '
          f'Isp = {ideal["Isp"]:.2f} s, c* = {ideal["cstar"]:.1f} m/s\n')

    solver = NozzleSolver(grid, gas, bcs, flux='roe', limiter='vanalbada',
                          viscous=not args.euler, cfl=args.cfl)
    solver.initialize(quasi1d.initial_field(grid, gas, p0, T0, p_amb))
    ckpt = out + '.npz'
    if os.path.exists(ckpt) and not args.fresh:
        solver.load(ckpt)
        print(f'resuming from {ckpt} at iteration {solver.iter}')
    solver.run(max_iter=args.iters, tol=args.tol, print_every=args.print_every,
               checkpoint_every=args.checkpoint_every, checkpoint_path=ckpt)

    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    solver.save(out + '.npz')
    post.report(solver, ideal)

    f = solver.fields()
    print(f'\nThroat Reynolds number : {post.throat_reynolds(solver):.0f}')
    print(f'Exit Mach, axis        : {f["M"][-1, 0]:.3f}')
    print(f'Boundary layer at exit : '
          f'{post.boundary_layer_thickness(solver)[-1] * 1e3:.4f} mm of '
          f'{contour.r_exit * 1e3:.4f} mm')
    print("\nHayn measured Isp ~63-75 s for cold N2 micronozzles at this scale "
          "(Fig. 73).")

    if not args.no_plots:
        from fvm import plots
        for p in plots.save_all(solver, prefix=out, ideal=ideal):
            print('  ' + p)


if __name__ == '__main__':
    main()
