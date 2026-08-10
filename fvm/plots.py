"""Matplotlib figures for the nozzle solution.

Colour choices follow one rule each: field contours use a single-hue ramp that
runs light to dark (never a rainbow -- perceptual non-uniformity in jet-style
maps invents features that are not in the data), and line plots draw from a
fixed, colour-vision-deficiency-checked categorical order rather than cycling.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# -- palette ---------------------------------------------------------------
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
         "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
         "#0d366b"]
_ORANGE = ["#fde4d7", "#fbcbb2", "#f8b18d", "#f49768", "#f07d4b", "#eb6834",
           "#d95926", "#bf4a1c", "#a13c15", "#82300f", "#64240a"]

CMAP_BLUE = LinearSegmentedColormap.from_list("seq_blue", _BLUE)
CMAP_ORANGE = LinearSegmentedColormap.from_list("seq_orange", _ORANGE)

#: Which one-hue ramp each field gets. Temperature takes the warm ramp; the
#: rest share the default blue.
FIELD_CMAP = {"T": CMAP_ORANGE}

FIELD_LABEL = {
    "M": "Mach number", "p": "Static pressure [bar]", "T": "Static temperature [K]",
    "rho": "Density [kg/m$^3$]", "u": "Axial velocity [m/s]",
    "v": "Radial velocity [m/s]", "speed": "Speed [m/s]",
}
FIELD_SCALE = {"p": 1e-5}


def apply_style():
    matplotlib.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK_2,
        "axes.titlecolor": INK, "axes.titlesize": 11, "axes.titleweight": "medium",
        "axes.labelsize": 9.5, "axes.linewidth": 0.8,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
        "axes.axisbelow": True,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": 8.5,
        "legend.labelcolor": INK_2,
        "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans", "sans-serif"],
        "figure.dpi": 130,
    })


def _mirror(ax, grid, C, cmap, norm=None, **kw):
    """Draw a cell field on both sides of the symmetry axis."""
    x, r = grid.x_n * 1e3, grid.r_n * 1e3
    m1 = ax.pcolormesh(x, r, C, cmap=cmap, norm=norm, shading="flat",
                       rasterized=True, **kw)
    ax.pcolormesh(x, -r, C, cmap=cmap, norm=norm, shading="flat",
                  rasterized=True, **kw)
    return m1


def _dress(ax, grid, title=None):
    ax.set_xlabel("Axial station, x [mm]")
    ax.set_ylabel("Radius, r [mm]")
    if title:
        ax.set_title(title, loc="left")
    ax.set_aspect("equal")
    ax.grid(False)
    rw = grid.r_n[:, -1] * 1e3
    ax.plot(grid.x_n[:, -1] * 1e3, rw, color=INK, lw=1.0)
    ax.plot(grid.x_n[:, -1] * 1e3, -rw, color=INK, lw=1.0)
    ax.axhline(0.0, color=AXIS, lw=0.6, ls=(0, (4, 4)))


# -- figures ---------------------------------------------------------------
def plot_grid(grid, stride=1, ax=None, title="Computational grid"):
    apply_style()
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 3.4))
    x, r = grid.x_n * 1e3, grid.r_n * 1e3
    for i in range(0, grid.ni + 1, stride):
        ax.plot(x[i], r[i], color=SERIES[0], lw=0.35, alpha=0.75)
        ax.plot(x[i], -r[i], color=SERIES[0], lw=0.35, alpha=0.75)
    for j in range(0, grid.nj + 1, stride):
        ax.plot(x[:, j], r[:, j], color=SERIES[0], lw=0.35, alpha=0.75)
        ax.plot(x[:, j], -r[:, j], color=SERIES[0], lw=0.35, alpha=0.75)
    _dress(ax, grid, f"{title} -- {grid.ni} x {grid.nj} cells")
    return ax.figure


def plot_field(solver, name="M", ax=None, log=False, xlim=None, focus=False,
               title=None):
    apply_style()
    g = solver.grid
    f = solver.fields()
    C = f[name] * FIELD_SCALE.get(name, 1.0)
    cmap = FIELD_CMAP.get(name, CMAP_BLUE)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 3.6))
    norm = matplotlib.colors.LogNorm(vmin=max(C.min(), C.max() * 1e-6),
                                     vmax=C.max()) if log else None
    m = _mirror(ax, g, C, cmap, norm=norm)
    cb = ax.figure.colorbar(m, ax=ax, pad=0.02, fraction=0.035)
    cb.set_label(FIELD_LABEL.get(name, name), color=INK_2, fontsize=9)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor(AXIS)
    _dress(ax, g, title or FIELD_LABEL.get(name, name))

    if xlim is None and focus:
        xlim = _focus_xlim(g)
    if xlim:
        ax.set_xlim(*xlim)
        sel = (g.x_n[:, -1] * 1e3 >= xlim[0]) & (g.x_n[:, -1] * 1e3 <= xlim[1])
        if sel.any():
            rmax = g.r_n[sel, -1].max() * 1e3 * 1.15
            ax.set_ylim(-rmax, rmax)
    return ax.figure


def _focus_xlim(grid):
    """Trim the constant-radius chamber, which carries no gradients."""
    c = grid.contour
    x0 = getattr(c, "x_conv_start", c.x_start)
    pad = 0.03 * (c.x_end - x0)
    return ((x0 - pad) * 1e3, (c.x_end + pad) * 1e3)


def plot_centerline(solver, ideal=None, axes=None, focus=True):
    """Axial distributions of Mach number and pressure against quasi-1-D."""
    from . import quasi1d, post
    apply_style()
    g = solver.grid
    f = solver.fields()
    x = g.xc[:, 0] * 1e3
    st = post.station_integrals(solver)

    if axes is None:
        _, axes = plt.subplots(2, 1, figsize=(8.4, 6.4), sharex=True)
    a0, a1 = axes

    ref = quasi1d.isentropic_field(g.xc[:, 0], g.contour, solver.gas,
                                   solver.bcs.p0, solver.bcs.T0)
    a0.plot(x, ref["M"], color=MUTED, lw=1.6, ls=(0, (5, 3)),
            label="Quasi-1-D isentropic")
    a0.plot(x, f["M"][:, 0], color=SERIES[0], label="CFD, axis")
    a0.plot(st["x"] * 1e3, st["M_avg"], color=SERIES[1],
            label="CFD, mass-averaged")
    a0.set_ylabel("Mach number")
    a0.legend(loc="upper left")
    a0.set_title("Mach number along the nozzle", loc="left")

    a1.semilogy(x, ref["p"] / 1e5, color=MUTED, lw=1.6, ls=(0, (5, 3)),
                label="Quasi-1-D isentropic")
    a1.semilogy(x, f["p"][:, 0] / 1e5, color=SERIES[0], label="CFD, axis")
    a1.semilogy(x, f["p"][:, -1] / 1e5, color=SERIES[1], label="CFD, wall")
    a1.axhline(solver.bcs.p_amb / 1e5, color=SERIES[3], lw=1.4, ls=(0, (1, 2)),
               label="Ambient")
    a1.set_ylabel("Static pressure [bar]")
    a1.set_xlabel("Axial station, x [mm]")
    a1.legend(loc="upper right")
    a1.set_title("Static pressure along the nozzle", loc="left")
    for a in (a0, a1):
        a.axvline(g.contour.x_throat * 1e3, color=AXIS, lw=0.8)
        if focus:
            a.set_xlim(*_focus_xlim(g))
    return a0.figure


def plot_wall(solver, axes=None, focus=True):
    """Wall pressure, skin friction, heat flux and grid resolution (y+)."""
    from . import post
    apply_style()
    w = post.wall_quantities(solver)
    x = w["x"] * 1e3

    if axes is None:
        _, axes = plt.subplots(2, 2, figsize=(9.6, 5.8), sharex=True)
    a = axes.ravel()

    a[0].semilogy(x, w["p"] / 1e5, color=SERIES[0])
    a[0].set_ylabel("Wall pressure [bar]")
    a[0].set_title("Wall static pressure", loc="left")

    a[1].plot(x, w["tau"], color=SERIES[1])
    a[1].set_ylabel(r"$\tau_w$ [Pa]")
    a[1].set_title("Wall shear stress", loc="left")

    # An adiabatic wall carries no heat flux by construction, so show the
    # recovery temperature the wall would have to survive instead.
    if solver.bcs.wall_is_adiabatic:
        a[2].plot(x, w["T"], color=SERIES[2])
        a[2].axhline(solver.bcs.T0, color=MUTED, lw=1.2, ls=(0, (4, 3)))
        a[2].annotate(f"$T_0$ = {solver.bcs.T0:.0f} K", xy=(x[2], solver.bcs.T0),
                      xytext=(0, 4), textcoords="offset points",
                      color=INK_2, fontsize=8)
        a[2].set_ylabel("Wall temperature [K]")
        a[2].set_title("Adiabatic wall (recovery) temperature", loc="left")
    else:
        a[2].plot(x, w["q"] * 1e-3, color=SERIES[2])
        a[2].set_ylabel("Heat flux [kW/m$^2$]")
        a[2].set_title("Wall heat flux (positive into wall)", loc="left")
    a[2].set_xlabel("Axial station, x [mm]")

    a[3].semilogy(x, np.maximum(w["y_plus"], 1e-6), color=SERIES[3])
    a[3].axhline(1.0, color=MUTED, lw=1.2, ls=(0, (4, 3)))
    a[3].annotate("y+ = 1", xy=(x[len(x) // 2], 1.0), xytext=(0, 4),
                  textcoords="offset points", color=INK_2, fontsize=8)
    a[3].set_ylabel("y$^+$")
    a[3].set_xlabel("Axial station, x [mm]")
    a[3].set_title("Wall-normal grid resolution", loc="left")

    for ax in a:
        ax.axvline(solver.grid.contour.x_throat * 1e3, color=AXIS, lw=0.8)
        if focus:
            ax.set_xlim(*_focus_xlim(solver.grid))
    return a[0].figure


def plot_exit_profile(solver, ax=None):
    """Radial profiles at the exit plane -- this is where the viscous loss lives."""
    apply_style()
    g = solver.grid
    f = solver.fields()
    r = g.rc[-1] * 1e3
    rw = g.r_n[-1, -1] * 1e3

    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(8.6, 4.0), sharey=True)
    ax[0].plot(f["M"][-1], r, color=SERIES[0], marker="o", ms=3.0)
    ax[0].set_xlabel("Mach number")
    ax[0].set_ylabel("Radius, r [mm]")
    ax[0].set_title("Exit-plane Mach profile", loc="left")

    # Show static and total temperature together. They differ by ~1600 K here,
    # and plotting only the static one invites the reading that the wall
    # carries more energy than the core -- it carries less, it just cannot
    # move, so none of its energy is kinetic.
    Tt = f["T"][-1] + (f["u"][-1] ** 2 + f["v"][-1] ** 2) / (2.0 * solver.gas.cp)
    ax[1].plot(f["T"][-1], r, color=SERIES[1], marker="o", ms=3.0,
               label="Static, $T$")
    ax[1].plot(Tt, r, color=SERIES[2], label="Total, $T_0 = T + q^2/2c_p$")
    ax[1].axvline(solver.bcs.T0, color=MUTED, lw=1.2, ls=(0, (4, 3)))
    ax[1].annotate(f"chamber $T_0$ = {solver.bcs.T0:.0f} K",
                   xy=(solver.bcs.T0, r[0]), xytext=(-6, 2),
                   textcoords="offset points", ha="right", rotation=90,
                   color=INK_2, fontsize=8)
    ax[1].set_xlabel("Temperature [K]")
    ax[1].set_title("Exit-plane temperature profile", loc="left")
    ax[1].legend(loc="lower left")
    for a in ax:
        a.axhline(rw, color=INK, lw=1.0)
        a.annotate("wall", xy=(a.get_xlim()[0], rw), xytext=(4, -10),
                   textcoords="offset points", color=INK_2, fontsize=8)
    return ax[0].figure


def plot_convergence(solver, ax=None):
    apply_style()
    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 3.6))
    h = np.array(solver.history)
    ax.semilogy(h[:, 0], h[:, 2], color=SERIES[0])
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Density residual, normalised")
    ax.set_title("Convergence history", loc="left")
    ax.annotate(f"final {h[-1, 2]:.2e}", xy=(h[-1, 0], h[-1, 2]),
                xytext=(-6, 8), textcoords="offset points",
                ha="right", color=INK_2, fontsize=8.5)
    return ax.figure


def plot_mass_conservation(solver, ax=None):
    """Mass flow through every axial station -- the cheapest convergence check."""
    from . import post
    apply_style()
    st = post.station_integrals(solver)
    if ax is None:
        _, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.plot(st["x"] * 1e3, st["mdot"] * 1e6, color=SERIES[0])
    ax.set_xlabel("Axial station, x [mm]")
    ax.set_ylabel("Mass flow [mg/s]")
    ax.set_title("Mass flow through each axial station", loc="left")
    ax.axvline(solver.grid.contour.x_throat * 1e3, color=AXIS, lw=0.8)
    return ax.figure


def save_all(solver, prefix="out/nozzle", ideal=None, dpi=140):
    """Write the standard figure set. Returns the list of paths written."""
    import os
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    paths = []
    figs = [
        ("grid", lambda: plot_grid(solver.grid)),
        ("mach", lambda: plot_field(solver, "M")),
        ("mach_nozzle", lambda: plot_field(
            solver, "M", focus=True, title="Mach number -- throat and divergent section")),
        ("pressure", lambda: plot_field(solver, "p", log=True)),
        ("temperature", lambda: plot_field(solver, "T", focus=True)),
        ("centerline", lambda: plot_centerline(solver, ideal)),
        ("wall", lambda: plot_wall(solver)),
        ("exit_profile", lambda: plot_exit_profile(solver)),
        ("mass_conservation", lambda: plot_mass_conservation(solver)),
    ]
    if solver.history:
        figs.append(("convergence", lambda: plot_convergence(solver)))
    for name, fn in figs:
        fig = fn()
        fig.tight_layout()
        path = f"{prefix}_{name}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths
