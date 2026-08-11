"""System-level coupling: catalyst bed to nozzle.

The loop
-------
Bed and nozzle are coupled through two facts that each depend on the other::

    chamber pressure = bed inlet pressure - bed pressure drop
    mass flow        = Cd * chamber pressure * A_throat / c*

The bed's pressure drop depends on mass flux, and ``c*`` depends on the bed
exit composition, which depends on the pressure the bed ran at. So neither can
be evaluated without the other and the system is solved by iteration.

It is a well-behaved fixed point. Pressure drop is a weak function of the
absolute pressure level -- the Ergun inertial term goes as ``G^2/rho`` and
density is proportional to pressure, so the drop is close to inversely
proportional to it -- which makes the residual nearly linear and a secant
method converge in a handful of steps.

Two directions, both physically meaningful:

* **Given mass flow** -- the design question. Size a thruster for a target
  flow and find what feed pressure it needs.
* **Given feed pressure** -- the operational question. A regulator sets the
  upstream pressure and the hardware decides the flow.

Cost, and where the expensive model goes
----------------------------------------
The iteration needs many bed integrations, so the nozzle side of the loop uses
the quasi-1-D isentropic relations rather than the CFD solver. The expensive
2-D viscous solve is run *once*, afterwards, on the converged chamber
conditions -- and its discharge coefficient can then be fed back through
``Cd`` to re-converge the system. That feedback is the point: at a 0.29 mm
throat ``Cd`` is around 0.93, and ignoring it puts the chamber pressure out by
the same 7%.
"""
import numpy as np

from . import quasi1d
from .plugflow import PlugFlowReactor
from .thermo import G0, PerfectGas


def perfect_gas_from_bed(exit_conditions, mechanism, Pr=0.7, omega=0.7):
    """Build a frozen-composition :class:`PerfectGas` from a bed exit state.

    The nozzle solver wants a single calorically perfect fluid, which is the
    right model downstream of the bed: residence time in the nozzle is far too
    short for the composition to shift, so it is frozen at whatever the bed
    delivered. Gamma is evaluated at the chamber temperature and held.
    """
    T = exit_conditions["T"]
    Y = exit_conditions["Y"]
    return PerfectGas(gamma=exit_conditions["gamma"], MW=exit_conditions["MW"],
                      Pr=Pr, mu_ref=float(mechanism.viscosity(T, Y)),
                      T_mu_ref=T, mu_law="power", omega=omega)


def vapor_region_inlet(mechanism, T_vapor, T_feed=298.15):
    """Composition where the vapour region begins, from energy conservation.

    The bed is fed liquid, and the energy to vaporise it comes from
    decomposing some of it. So the extent of that pre-decomposition is not a
    free parameter -- it is fixed by requiring

        h(T_vapor, Y) = h_liquid(T_feed)

    Getting this wrong is not subtle. Starting the integration from *pure*
    hydrazine vapour at ``T_vapor`` hands the bed the vaporisation enthalpy
    for nothing, and the adiabatic flame temperature comes out hundreds of
    kelvin too high -- 1428 K at complete dissociation against the correct
    868 K.

    Returns the mass fractions after a fraction ``f`` of the hydrazine has
    decomposed by ``N2H4 -> 4/3 NH3 + 1/3 N2``, with ``f`` solved for.
    """
    mix = mechanism.mixture
    h_target = mechanism.inlet_enthalpy(T_feed)

    def composition(f):
        return mix.from_moles({"N2H4": 1.0 - f,
                               "NH3": 4.0 / 3.0 * f,
                               "N2": f / 3.0,
                               "H2": 0.0})

    def residual(f):
        return float(mix.h(T_vapor, composition(f))) - h_target

    # residual decreases with f: decomposition is exothermic, so at a fixed
    # temperature a more-decomposed mixture carries less enthalpy.
    lo, hi = 0.0, 1.0
    if residual(lo) < 0.0:
        raise ValueError(
            f"even undecomposed vapour at T_vapor = {T_vapor:.1f} K carries "
            f"less enthalpy than the liquid feed; T_vapor is too low")
    if residual(hi) > 0.0:
        raise ValueError(
            f"decomposing all the hydrazine still leaves more enthalpy than "
            f"the feed at T_vapor = {T_vapor:.1f} K; T_vapor is too high")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    f = 0.5 * (lo + hi)
    return composition(f), f


class ThrusterSolution:
    """A converged bed + nozzle operating point."""

    def __init__(self, bed_solution, gas, contour, p_ambient, mdot, Cd, ideal):
        self.bed = bed_solution
        self.gas = gas
        self.contour = contour
        self.p_ambient = p_ambient
        self.mdot = mdot
        self.Cd = Cd
        self.ideal = ideal
        self.chamber = bed_solution.exit_conditions()

    @property
    def p_feed(self):
        return float(self.bed.p[0])

    @property
    def p_chamber(self):
        return float(self.bed.p[-1])

    @property
    def bed_pressure_drop(self):
        return self.p_feed - self.p_chamber

    @property
    def thrust(self):
        return self.ideal["thrust"]

    @property
    def Isp(self):
        return self.ideal["Isp"]

    def report(self, quiet=False):
        c, b = self.chamber, self.bed
        L = ["", "Thruster System -- bed and nozzle coupled", "-" * 66]
        L.append(f"{'Feed pressure':<32}{self.p_feed / 1e5:>12.3f}  bar")
        L.append(f"{'Bed pressure drop':<32}{self.bed_pressure_drop / 1e5:>12.3f}  bar"
                 f"   ({100 * self.bed_pressure_drop / self.p_feed:.1f}% of feed)")
        L.append(f"{'Chamber pressure':<32}{self.p_chamber / 1e5:>12.3f}  bar")
        L.append("-" * 66)
        L.append(f"{'Bed exit temperature':<32}{c['T']:>12.1f}  K")
        L.append(f"{'Ammonia dissociation X':<32}{c['dissociation']:>12.3f}")
        L.append(f"{'Exit molecular weight':<32}{c['MW']:>12.3f}  kg/kmol")
        L.append(f"{'Exit gamma':<32}{c['gamma']:>12.4f}")
        L.append(f"{'Characteristic velocity':<32}{c['cstar']:>12.1f}  m/s")
        L.append(f"{'Catalyst peak temperature':<32}{b.T_solid.max():>12.1f}  K")
        L.append("-" * 66)
        L.append(f"{'Mass flow':<32}{self.mdot * 1e3:>12.5f}  g/s")
        L.append(f"{'Discharge coefficient':<32}{self.Cd:>12.4f}")
        L.append(f"{'Thrust':<32}{self.thrust * 1e3:>12.3f}  mN")
        L.append(f"{'Specific impulse':<32}{self.Isp:>12.2f}  s")
        L.append(f"{'Exit Mach':<32}{self.ideal['M_exit']:>12.3f}")
        L.append("-" * 66)
        text = "\n".join(L)
        if not quiet:
            print(text)
        return text


class ThrusterSystem:
    """Couples a catalyst bed to a nozzle and solves for the operating point.

    Parameters
    ----------
    mechanism, bed, contour
    p_ambient : float
        Back pressure [Pa].
    T_vapor : float
        Temperature at which the vapour region begins [K]. See
        :mod:`fvm.plugflow` -- vaporisation is not modelled, so this is a
        modelling input rather than a result, and the answer is sensitive to
        it through the Arrhenius rates.
    Cd : float
        Throat discharge coefficient. 1.0 is the inviscid value; run the CFD
        solver and feed the result back for a real thruster.
    """

    def __init__(self, mechanism, bed, contour, p_ambient=0.0,
                 T_vapor=455.6, T_feed=298.15, Cd=1.0, Y_inlet=None,
                 reactor_kw=None):
        self.mechanism = mechanism
        self.bed = bed
        self.contour = contour
        self.p_ambient = float(p_ambient)
        self.T_vapor = float(T_vapor)
        self.T_feed = float(T_feed)
        self.Cd = float(Cd)
        if Y_inlet is None:
            Y_inlet, self.pre_decomposition = vapor_region_inlet(
                mechanism, self.T_vapor, T_feed)
        else:
            self.pre_decomposition = float("nan")
        self.Y_inlet = Y_inlet
        self.reactor = PlugFlowReactor(mechanism, bed, **(reactor_kw or {}))

    @property
    def throat_area(self):
        return np.pi * self.contour.r_throat ** 2

    # -- one pass through the physics -------------------------------------
    def _evaluate(self, p_inlet, G, n_output=3, profiles=False, rtol=1e-6):
        """Integrate the bed and work out what mass flow the throat would pass.

        Defaults to exit-state-only: the iteration needs one number from each
        pass, and reconstructing full profiles costs a bracketed catalyst
        temperature solve at every station.
        """
        sol = self.reactor.solve(G=G, p_inlet=p_inlet, T_inlet=self.T_vapor,
                                 Y_inlet=self.Y_inlet, n_output=n_output,
                                 profiles=profiles, rtol=rtol)
        exit_c = sol.exit_conditions()
        mdot_throat = (self.Cd * exit_c["p"] * self.throat_area / exit_c["cstar"])
        return sol, exit_c, mdot_throat

    # -- design case: mass flow given -------------------------------------
    def solve_for_mdot(self, mdot, tol=1e-6, max_iter=40, n_output=200,
                       rtol=1e-8):
        """Find the feed pressure that passes ``mdot`` through bed and throat.

        Solved as a fixed point rather than by bracketing::

            p_feed <- mdot c* / (Cd A_t) + dP_bed(p_feed)

        The first term is exact and the correction is weak, so this converges
        in a handful of bed integrations where bisection needed forty. Each
        integration is the expensive part, so iteration count is the whole
        cost story.

        ``tol`` cannot be pushed below the ODE integrator's own noise floor:
        asking for 1e-8 on the pressure while each bed pass carries 1e-5
        relative error simply never converges. ``rtol`` is therefore set
        tighter than the outer tolerance rather than looser.
        """
        G = self.bed.mass_flux(mdot)
        _, exit_c, _ = self._evaluate(2e5, G, rtol=rtol)
        p = mdot * exit_c["cstar"] / (self.Cd * self.throat_area)

        for _ in range(max_iter):
            sol, exit_c, _ = self._evaluate(p, G, rtol=rtol)
            drop = float(sol.p[0] - sol.p[-1])
            p_new = mdot * exit_c["cstar"] / (self.Cd * self.throat_area) + drop
            converged = abs(p_new - p) < tol * p_new
            p = p_new
            if converged:
                break
        else:
            raise RuntimeError("feed pressure iteration did not converge")

        sol, exit_c, _ = self._evaluate(p, G, n_output=n_output, profiles=True)
        return self._finish(sol, exit_c, mdot)

    # -- operational case: feed pressure given ----------------------------
    def solve_for_feed_pressure(self, p_feed, tol=1e-6, max_iter=60,
                                n_output=200, rtol=1e-8):
        """Find the mass flow a given feed pressure establishes.

        Fixed point on mass flow: the throat sets what the bed exit pressure
        can pass, and the bed drop sets that exit pressure.
        """
        _, exit_c, _ = self._evaluate(p_feed, self.bed.mass_flux(1e-6), rtol=rtol)
        mdot = self.Cd * p_feed * self.throat_area / exit_c["cstar"]

        for _ in range(max_iter):
            G = self.bed.mass_flux(mdot)
            sol, exit_c, mdot_throat = self._evaluate(p_feed, G, rtol=rtol)
            # Under-relax: the bed drop responds to mdot roughly quadratically,
            # so a full step overshoots and can oscillate.
            mdot_new = 0.5 * (mdot + mdot_throat)
            converged = abs(mdot_new - mdot) < tol * mdot_new
            mdot = mdot_new
            if converged:
                break
        else:
            raise RuntimeError("mass flow iteration did not converge")

        sol, exit_c, _ = self._evaluate(p_feed, self.bed.mass_flux(mdot),
                                        n_output=n_output, profiles=True)
        return self._finish(sol, exit_c, mdot)

    # -- shared tail ------------------------------------------------------
    def _finish(self, sol, exit_c, mdot):
        gas = perfect_gas_from_bed(exit_c, self.mechanism)
        ideal = quasi1d.ideal_performance(self.contour, gas, exit_c["p"],
                                          exit_c["T"], self.p_ambient)
        # Replace the isentropic mass flow with the converged system value, so
        # thrust and Isp reflect the actual operating point rather than the
        # ideal choked flow at this chamber pressure.
        ue = ideal["u_exit"]
        Ae = ideal["exit_area"]
        F = mdot * ue + (ideal["p_exit"] - self.p_ambient) * Ae
        ideal = dict(ideal)
        ideal.update(mdot=mdot, thrust=F, Isp=F / (mdot * G0),
                     Cf=F / (exit_c["p"] * self.throat_area))
        return ThrusterSolution(sol, gas, self.contour, self.p_ambient,
                                mdot, self.Cd, ideal)
