"""Solve the complete thruster: liquid hydrazine -> catalyst bed -> nozzle.

Couples the reacting plug-flow bed model to the nozzle, iterating until the
chamber pressure, the bed pressure drop and the choked throat are mutually
consistent. Edit the input block below, or override from the command line::

    python run_thruster.py --mdot 0.0443          # design: flow given
    python run_thruster.py --p-feed 8.5           # operation: pressure given
    python run_thruster.py --cd 0.93              # feed back a CFD discharge coefficient
    python run_thruster.py --n-h2 1.6 --a-nh3 0.3e11   # alternative kinetics

The nozzle side of the iteration uses the quasi-1-D isentropic relations,
because the loop needs many bed integrations. Run ``run_fvm_nozzle.py`` once
afterwards on the chamber conditions printed here to get the real discharge
coefficient and thrust efficiency, then pass that Cd back in with ``--cd``.
"""
import argparse
import os

import numpy as np
import metronos as mt

from fvm import (CatalystBed, ConicalNozzle, HydrazineShell405,
                 PackedSpheres, ReticulatedFoam)
from fvm.thruster import ThrusterSystem, vapor_region_inlet

# ── Inputs ───────────────────────────────────────────────────────────────────
d_c   = mt.Quantity(8.3,  'mm')     # Chamber / bed diameter
d_t   = mt.Quantity(0.29, 'mm')     # Throat diameter
expan = 100.0                       # Nozzle area ratio
L_bed = mt.Quantity(20.0, 'mm')     # Catalyst bed length
P_a   = mt.Quantity(5.0,  'Torr')   # Ambient pressure

T_feed  = mt.Quantity(298.15, 'K')  # Liquid hydrazine feed temperature
T_vapor = mt.Quantity(455.6,  'K')  # Where the vapour region begins

theta_conv, theta_div = 35.0, 15.0


def build_bed(args):
    d = d_c.to('m').value
    L = (args.bed_length * 1e-3) if args.bed_length else L_bed.to('m').value
    if args.packing == 'foam':
        return CatalystBed.uniform(d, L, ReticulatedFoam(ppi=args.ppi))
    return CatalystBed.uniform(
        d, L, PackedSpheres(mesh=(args.mesh_coarse, args.mesh_fine),
                            eps=args.voidage))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    op = ap.add_mutually_exclusive_group()
    op.add_argument('--mdot', type=float, default=None, metavar='G_S',
                    help='design case: solve for the feed pressure this flow needs [g/s]')
    op.add_argument('--p-feed', type=float, default=None, metavar='BAR',
                    help='operational case: solve for the flow this feed pressure gives')

    ap.add_argument('--cd', type=float, default=1.0,
                    help='throat discharge coefficient; feed back from run_fvm_nozzle.py')
    ap.add_argument('--packing', default='spheres', choices=['spheres', 'foam'])
    ap.add_argument('--mesh-coarse', type=int, default=25)
    ap.add_argument('--mesh-fine', type=int, default=30)
    ap.add_argument('--voidage', type=float, default=0.38)
    ap.add_argument('--ppi', type=float, default=60.0)
    ap.add_argument('--bed-length', type=float, default=None, metavar='MM')
    ap.add_argument('--t-vapor', type=float, default=None, metavar='K')

    ap.add_argument('--n-h2', type=float, default=1.0,
                    help='hydrogen inhibition order; 1.0 reproduces Kesten, '
                         '1.6 is his published Eq. 43')
    ap.add_argument('--a-nh3', type=float, default=1.0e11,
                    help='ammonia pre-exponential in Kesten units')
    ap.add_argument('--out', default='out/thruster')
    ap.add_argument('--no-plots', action='store_true')
    args = ap.parse_args()

    if args.mdot is None and args.p_feed is None:
        args.mdot = 0.0443              # g/s, roughly 100 mN

    np.seterr(all='ignore')
    T_vap = args.t_vapor if args.t_vapor else T_vapor.to('K').value

    mech = HydrazineShell405(nH2=args.n_h2, A_NH3=args.a_nh3)
    bed = build_bed(args)
    nozzle = ConicalNozzle(d_c.to('m').value / 2, d_t.to('m').value / 2,
                           expan, 4e-3, theta_conv=theta_conv, theta_div=theta_div)
    system = ThrusterSystem(mech, bed, nozzle, p_ambient=P_a.to('Pa').value,
                            T_vapor=T_vap, T_feed=T_feed.to('K').value, Cd=args.cd)

    print(bed.summary())
    _, f = vapor_region_inlet(mech, T_vap, T_feed.to('K').value)
    print(f"\nVapour region begins at {T_vap:.1f} K with {100 * f:.1f}% of the "
          f"hydrazine already decomposed.")
    print("That extent is not assumed -- it is what energy conservation demands "
          "to vaporise the feed.")
    print(f"Kinetics: hydrogen order n = {args.n_h2}, A_NH3 = {args.a_nh3:.2e} "
          f"(Kesten units), Cd = {args.cd}")

    if args.mdot is not None:
        print(f"\nSolving for the feed pressure that passes {args.mdot:g} g/s ...")
        sol = system.solve_for_mdot(args.mdot * 1e-3)
    else:
        print(f"\nSolving for the flow established by {args.p_feed:g} bar ...")
        sol = system.solve_for_feed_pressure(args.p_feed * 1e5)

    sol.report()

    loading = bed.loading_lb_in2_s(sol.mdot)
    print(f"{'Bed loading':<32}{loading:>12.4f}  lb/in^2-s")
    if loading < 0.02:
        print("  Below the 0.02-0.05 lb/in^2-s band typical of Shell 405 practice:")
        print("  conservative for decomposition, but the bed is larger than it needs to be.")
    elif loading > 0.05:
        print("  Above typical Shell 405 practice -- check for incomplete decomposition.")
    print("-" * 66)
    print("\nHand these chamber conditions to run_fvm_nozzle.py for the real")
    print("discharge coefficient, then re-run with --cd to close the loop:")
    c = sol.chamber
    print(f"  p0 = {c['p'] / 1e5:.4f} bar   T0 = {c['T']:.1f} K   "
          f"MW = {c['MW']:.3f}   gamma = {c['gamma']:.4f}")

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out + '_report.txt', 'w', encoding='utf-8') as fh:
        fh.write(sol.report(quiet=True))
    print(f"\nreport written to {args.out}_report.txt")

    if not args.no_plots:
        from fvm import plots
        path = plots.plot_bed_profiles(sol.bed, prefix=args.out)
        print(f"bed profiles written to {path}")


if __name__ == '__main__':
    main()
