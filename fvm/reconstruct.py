"""MUSCL reconstruction and slope limiters.

Reconstruction is performed on primitive variables ``[rho, u, v, p]``, which
keeps the limiter monotone in pressure and density and makes the positivity
fallback trivial.
"""
import numpy as np

EPS = 1e-30


def lim_none(a, b):
    return np.zeros_like(a)


def lim_minmod(a, b):
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))


def lim_vanleer(a, b):
    s = a * b
    return np.where(s > 0.0, 2.0 * s / np.where(np.abs(a + b) < EPS, EPS, a + b), 0.0)


def lim_vanalbada(a, b):
    """van Albada, switched off across an extremum.

    The bare averaged-slope expression is *not* TVD: for a = -2, b = 3 it
    returns -0.46 rather than 0, so it reconstructs a slope straight through a
    local extremum. The ``a*b > 0`` gate is what makes it a limiter.
    """
    e = 1e-12
    s = (a * (b * b + e) + b * (a * a + e)) / (a * a + b * b + 2.0 * e)
    return np.where(a * b > 0.0, s, 0.0)


def lim_mc(a, b):
    """Monotonized central: sharpest of the classic second-order limiters."""
    return (0.5 * (np.sign(a) + np.sign(b))
            * np.minimum(np.minimum(2.0 * np.abs(a), 2.0 * np.abs(b)),
                         0.5 * np.abs(a + b)))


LIMITERS = {
    "none": lim_none,
    "minmod": lim_minmod,
    "vanleer": lim_vanleer,
    "vanalbada": lim_vanalbada,
    "mc": lim_mc,
}


def get_limiter(name):
    try:
        return LIMITERS[name.lower()]
    except KeyError:
        raise ValueError(f"unknown limiter {name!r}; "
                         f"choose from {sorted(LIMITERS)}") from None


def _sl(ndim, axis, a, b):
    return tuple(slice(None) if k != axis else slice(a, b) for k in range(ndim))


def muscl(W, axis, ng, n, limiter):
    """Reconstruct left/right face states along ``axis``.

    ``W`` has ``ng`` ghost layers on each side of ``n`` interior cells along
    ``axis``. Returns ``(WL, WR)`` for the ``n+1`` faces of that family.

    Faces where the reconstruction would produce a non-positive density or
    pressure fall back to first order locally, which is what keeps the solver
    alive through the throat and across startup transients.
    """
    d = W.ndim
    im2 = W[_sl(d, axis, ng - 2, ng + n - 1)]
    im1 = W[_sl(d, axis, ng - 1, ng + n)]
    ip0 = W[_sl(d, axis, ng, ng + n + 1)]
    ip1 = W[_sl(d, axis, ng + 1, ng + n + 2)]

    WL = im1 + 0.5 * limiter(im1 - im2, ip0 - im1)
    WR = ip0 - 0.5 * limiter(ip0 - im1, ip1 - ip0)

    badL = (WL[0] <= 0.0) | (WL[3] <= 0.0)
    badR = (WR[0] <= 0.0) | (WR[3] <= 0.0)
    WL = np.where(badL[None, ...], im1, WL)
    WR = np.where(badR[None, ...], ip0, WR)
    return WL, WR
