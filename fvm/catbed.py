"""Catalyst bed geometry, packing and pressure drop.

What the bed model has to supply
--------------------------------
The plug-flow solver needs three things at every axial station: how much
catalyst surface is available per unit volume (``a_v``), how big the particles
are (``d_p``, which sets the mass-transfer film thickness), and how fast
pressure is falling. Everything here exists to provide those.

Beds are built from *sections*, because real ones are graded. Kesten's standard
configuration is 25-30 mesh granules for the first 0.2 in and 1/8 x 1/8 in
cylindrical pellets for the remainder: fine material where the reaction is
fastest and surface area matters, coarse material downstream where pressure
drop would otherwise dominate. A single-section bed is just the degenerate
case.

Pressure drop
-------------
Both packed beds and foams use the Ergun form,

    -dp/dz = 150 mu (1-eps)^2 / (eps^3 d^2) u  +  1.75 rho (1-eps) / (eps^3 d) u^2

with ``u`` the superficial velocity ``G/rho`` and ``d`` the effective particle
diameter ``sphericity * d_volume_equivalent``. The first term is viscous
(Darcy, linear in u), the second inertial (Forchheimer, quadratic). In a
hydrazine bed the gas accelerates by orders of magnitude as it decomposes and
heats, so the inertial term ends up dominating even though bed Reynolds
numbers start low -- which is why the quadratic term cannot be dropped.

Foams are the same equation with a much higher void fraction and a pore-scale
length, which is standard practice; a genuinely separate Forchheimer
correlation with a fitted permeability would be better if you have foam
pressure-drop data to fit it to.
"""
import numpy as np

#: US standard sieve openings [m], for the mesh sizes catalyst is sold in.
SIEVE_OPENING = {
    8: 2.36e-3, 10: 2.00e-3, 12: 1.70e-3, 14: 1.40e-3, 16: 1.18e-3,
    18: 1.00e-3, 20: 850e-6, 25: 710e-6, 30: 600e-6, 35: 500e-6,
    40: 425e-6, 45: 355e-6, 50: 300e-6, 60: 250e-6, 70: 212e-6,
    80: 180e-6, 100: 150e-6,
}


def mesh_to_diameter(coarse, fine=None):
    """Particle diameter [m] for a sieve cut.

    ``mesh_to_diameter(25, 30)`` is material passing a 25 mesh and retained on
    a 30 mesh, i.e. between 0.600 and 0.710 mm; the geometric mean is returned
    as the representative size.
    """
    try:
        d_hi = SIEVE_OPENING[coarse]
    except KeyError:
        raise KeyError(f"no sieve opening for {coarse} mesh; "
                       f"have {sorted(SIEVE_OPENING)}") from None
    if fine is None:
        return d_hi
    try:
        d_lo = SIEVE_OPENING[fine]
    except KeyError:
        raise KeyError(f"no sieve opening for {fine} mesh") from None
    if d_lo >= d_hi:
        raise ValueError("the finer mesh number must be the larger number")
    return float(np.sqrt(d_hi * d_lo))


# ---------------------------------------------------------------------------
# packing media
# ---------------------------------------------------------------------------
class BedMedium:
    """Base class: a uniform packing with a void fraction and a length scale.

    Subclasses set ``eps`` (void fraction), ``d_v`` (diameter of the sphere of
    equal volume) and ``sphericity``. Everything else follows.
    """

    eps = 0.4
    d_v = 1e-3
    sphericity = 1.0
    name = "medium"

    @property
    def d_p(self):
        """Effective particle diameter [m] -- what the correlations want."""
        return self.sphericity * self.d_v

    @property
    def a_v(self):
        """Catalyst surface per unit *bed* volume [1/m].

        For a packing of convex particles, ``a_v = 6 (1 - eps) / d_p``. This is
        external geometric area only; the internal area of a porous support is
        far larger and is what a BET measurement reports, but film diffusion
        sees the external surface, which is what the mass-transfer rate needs.
        """
        return 6.0 * (1.0 - self.eps) / self.d_p

    def pressure_gradient(self, G, rho, mu):
        """``-dp/dz`` [Pa/m], positive for flow in +z. Ergun."""
        u = G / np.maximum(rho, 1e-12)
        e3 = self.eps ** 3
        viscous = 150.0 * mu * (1.0 - self.eps) ** 2 / (e3 * self.d_p ** 2) * u
        inertial = 1.75 * rho * (1.0 - self.eps) / (e3 * self.d_p) * u * u
        return viscous + inertial

    def reynolds(self, G, mu):
        """Particle Reynolds number, ``G d_p / (mu (1 - eps))``.

        The ``(1 - eps)`` normalisation is the packed-bed convention: below
        about 10 the viscous term dominates, above a few hundred the inertial
        one does.
        """
        return G * self.d_p / (np.maximum(mu, 1e-12) * (1.0 - self.eps))

    def __repr__(self):
        return (f"{type(self).__name__}(d_p={self.d_p * 1e3:.3f} mm, "
                f"eps={self.eps:.3f}, a_v={self.a_v:.0f} 1/m)")


class PackedSpheres(BedMedium):
    """Spherical or near-spherical granules."""

    def __init__(self, diameter=None, mesh=None, eps=0.38, sphericity=1.0,
                 name="spheres"):
        if (diameter is None) == (mesh is None):
            raise ValueError("give exactly one of diameter or mesh")
        if mesh is not None:
            diameter = mesh_to_diameter(*np.atleast_1d(mesh))
        self.d_v = float(diameter)
        self.eps = float(eps)
        self.sphericity = float(sphericity)
        self.name = name


class PackedCylinders(BedMedium):
    """Cylindrical pellets, e.g. Kesten's 1/8 x 1/8 in.

    Sphericity is computed from the geometry rather than tabulated: it is the
    ratio of the equal-volume sphere's area to the cylinder's own area, which
    for L = D works out at 0.874.
    """

    def __init__(self, diameter, length=None, eps=0.4, name="cylinders"):
        d, L = float(diameter), float(length if length is not None else diameter)
        volume = np.pi * d * d * L / 4.0
        area = np.pi * d * d / 2.0 + np.pi * d * L
        self.d_v = (6.0 * volume / np.pi) ** (1.0 / 3.0)
        self.sphericity = float(np.pi * self.d_v ** 2 / area)
        self.eps = float(eps)
        self.diameter, self.length = d, L
        self.name = name


class ReticulatedFoam(BedMedium):
    """Open-cell foam, characterised by pores per inch.

    Specific surface area uses the correlation already in
    ``thruster_sizing.py``::

        a_v = 1300 (PPI/20)^1.3     [m^2/m^3]

    which is an empirical fit, not a geometric identity -- unlike the packed
    bed, ``a_v`` here is *not* derived from ``eps`` and ``d_p``, so the two can
    disagree. Fit the correlation to your own foam if the surface area matters.
    """

    def __init__(self, ppi, eps=0.90, a_v_coefficient=1300.0,
                 a_v_exponent=1.3, name="foam"):
        self.ppi = float(ppi)
        self.eps = float(eps)
        self._a_v = a_v_coefficient * (self.ppi / 20.0) ** a_v_exponent
        # Pore diameter from pore count; the strut-scale length that sets the
        # pressure drop and the mass-transfer film.
        self.d_v = 0.0254 / self.ppi
        self.sphericity = 1.0
        self.name = name

    @property
    def a_v(self):
        return self._a_v


# ---------------------------------------------------------------------------
# the bed
# ---------------------------------------------------------------------------
class CatalystBed:
    """A cylindrical bed built from one or more axial sections.

    Parameters
    ----------
    diameter : float
        Bed diameter [m].
    sections : sequence of (length, medium)
        Ordered from the injector face downstream.
    """

    def __init__(self, diameter, sections):
        self.diameter = float(diameter)
        if not sections:
            raise ValueError("a bed needs at least one section")
        self.sections = [(float(L), m) for L, m in sections]
        if any(L <= 0 for L, _ in self.sections):
            raise ValueError("section lengths must be positive")
        self._edges = np.concatenate([[0.0], np.cumsum([L for L, _ in self.sections])])

    # -- geometry ---------------------------------------------------------
    @property
    def length(self):
        return float(self._edges[-1])

    @property
    def area(self):
        """Bed cross-sectional area [m^2]."""
        return 0.25 * np.pi * self.diameter ** 2

    @property
    def volume(self):
        return self.area * self.length

    def medium_at(self, z):
        """The packing at axial station ``z`` [m]."""
        i = int(np.clip(np.searchsorted(self._edges, z, side="right") - 1,
                        0, len(self.sections) - 1))
        return self.sections[i][1]

    def profile(self, attr, z):
        """Sample a medium property along z, e.g. ``profile('a_v', z)``."""
        z = np.atleast_1d(z)
        return np.array([getattr(self.medium_at(zi), attr) for zi in z])

    @property
    def catalyst_area(self):
        """Total external catalyst surface in the bed [m^2]."""
        return float(sum(L * m.a_v for L, m in self.sections)) * self.area

    # -- operating point --------------------------------------------------
    def mass_flux(self, mdot):
        """Superficial mass flux G [kg/(m^2 s)] -- the bed loading."""
        return mdot / self.area

    def loading_lb_in2_s(self, mdot):
        """Bed loading in the units hydrazine engines are usually specified in."""
        return self.mass_flux(mdot) * 0.00142233

    def residence_time(self, mdot, rho):
        """Crude residence time [s] using the void volume and a mean density.

        Only as good as the density estimate, which varies by orders of
        magnitude across a decomposing bed -- use the plug-flow solution if the
        number matters.
        """
        eps_mean = sum(L * m.eps for L, m in self.sections) / self.length
        return eps_mean * self.volume * rho / mdot

    def pressure_gradient(self, z, G, rho, mu):
        return self.medium_at(z).pressure_gradient(G, rho, mu)

    def pressure_drop(self, G, rho, mu, n=400):
        """Integrate Ergun over the bed for constant properties [Pa].

        ``rho`` and ``mu`` may be callables of z, which is how to get a
        meaningful answer for a reacting bed: density falls by orders of
        magnitude as liquid becomes hot gas, and the inertial term goes as
        1/rho, so assuming a constant density underestimates badly.
        """
        z = np.linspace(0.0, self.length, n)
        rho_f = rho if callable(rho) else (lambda _z: rho)
        mu_f = mu if callable(mu) else (lambda _z: mu)
        grad = np.array([self.pressure_gradient(zi, G, rho_f(zi), mu_f(zi))
                         for zi in z])
        return float(np.trapezoid(grad, z))

    # -- named configurations ---------------------------------------------
    @classmethod
    def kesten_standard(cls, diameter, length=0.25 * 0.3048):
        """Kesten's "standard bed configuration" (F910461-12, p.15).

        25-30 mesh granules for the first 0.2 in, 1/8 x 1/8 in cylindrical
        pellets for the remainder.
        """
        fine_length = 0.2 * 0.0254
        if length <= fine_length:
            raise ValueError("bed shorter than the 0.2 in fine-particle entry zone")
        return cls(diameter, [
            (fine_length, PackedSpheres(mesh=(25, 30), eps=0.38, name="25-30 mesh")),
            (length - fine_length,
             PackedCylinders(0.125 * 0.0254, 0.125 * 0.0254, eps=0.40,
                             name="1/8 x 1/8 in pellets")),
        ])

    @classmethod
    def uniform(cls, diameter, length, medium):
        return cls(diameter, [(length, medium)])

    def summary(self):
        L = [f"CatalystBed  D = {self.diameter * 1e3:.2f} mm, "
             f"L = {self.length * 1e3:.2f} mm, "
             f"A = {self.area * 1e6:.2f} mm^2"]
        z0 = 0.0
        for length, m in self.sections:
            L.append(f"  {z0 * 1e3:7.2f} - {(z0 + length) * 1e3:7.2f} mm : "
                     f"{m.name:<22s} d_p={m.d_p * 1e6:7.1f} um  "
                     f"eps={m.eps:.3f}  a_v={m.a_v:8.0f} 1/m")
            z0 += length
        L.append(f"  total external catalyst area: {self.catalyst_area * 1e4:.2f} cm^2")
        return "\n".join(L)

    def __repr__(self):
        return (f"CatalystBed(D={self.diameter * 1e3:.2f} mm, "
                f"L={self.length * 1e3:.2f} mm, {len(self.sections)} sections)")
