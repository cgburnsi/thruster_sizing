"""Interactive tool for fitting catalyst bed kinetics to measured data.

Run it with no arguments to load Kesten's reference case as a worked example::

    python fit_kinetics.py
    python fit_kinetics.py --data mybed.csv --bed-length 20 --bed-dia 8.3 --mdot 0.0443

The point of the tool is not to produce a number. It is to show whether your
data can support one. Two things it deliberately puts in front of you:

**The fitting window.** This model starts at the vapour-region inlet and does
not represent the liquid or liquid-vapour zones. Measurements from there cannot
be matched by any parameter value, and including them simply drags the fit off.
On the built-in Kesten case, excluding the first millimetre takes the
temperature residual from 59 K to 1.5 K.

**Sensitivity.** ``A_NH3`` and ``n_H2`` both control ammonia dissociation and
trade off against each other, so sparse data pins their combination while
leaving each poorly determined. The sensitivity readout says how much the
objective actually moves for a 10% change in each; a small number means that
parameter is unconstrained whatever value the fit returns.
"""
import argparse
import queue
import threading
import tkinter as tk
from tkinter import filedialog, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from fvm import plots
from fvm.catbed import CatalystBed, PackedSpheres
from fvm.fitting import PARAMETERS, BedData, FitCase, fit

SLIDER_ORDER = ("A_NH3", "n_H2", "T_vapor")


class FitApp:
    def __init__(self, root, case):
        self.root = root
        self.case = case
        self.queue = queue.Queue()
        self.busy = False
        self.pending = False
        self.solution = None
        root.title("Catalyst bed kinetics — fitting")

        self.values = {n: tk.DoubleVar(value=self._to_slider(n, PARAMETERS[n]["default"]))
                       for n in SLIDER_ORDER}
        self.readout = {n: tk.StringVar() for n in SLIDER_ORDER}
        self.status = tk.StringVar(value="ready")
        self.metrics = tk.StringVar(value="")
        self.zmin = tk.DoubleVar(value=1.0)
        self.zmax = tk.DoubleVar(value=case.data.z[-1] * 1e3)

        self._build(root)
        self._sync_readouts()
        self.request_solve()
        root.after(60, self._poll)

    # -- parameter <-> slider mapping (log where appropriate) -------------
    @staticmethod
    def _to_slider(name, value):
        return np.log10(value) if PARAMETERS[name]["log"] else value

    @staticmethod
    def _from_slider(name, value):
        return 10.0 ** value if PARAMETERS[name]["log"] else value

    def params(self):
        return {n: self._from_slider(n, self.values[n].get()) for n in SLIDER_ORDER}

    # -- layout -----------------------------------------------------------
    def _build(self, root):
        left = ttk.Frame(root, padding=8)
        left.grid(row=0, column=0, sticky="ns")
        right = ttk.Frame(root, padding=4)
        right.grid(row=0, column=1, sticky="nsew")
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)

        ttk.Label(left, text=f"Data: {self.case.data.name}",
                  font=("", 10, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
        ttk.Button(left, text="Load CSV…", command=self.load_csv).grid(
            row=1, column=0, columnspan=2, sticky="ew", pady=(2, 8))

        ttk.Label(left, text="Fitting window [mm]").grid(row=2, column=0,
                                                         columnspan=2, sticky="w")
        wf = ttk.Frame(left)
        wf.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        for i, (lbl, var) in enumerate((("from", self.zmin), ("to", self.zmax))):
            ttk.Label(wf, text=lbl).grid(row=0, column=2 * i)
            e = ttk.Entry(wf, textvariable=var, width=8)
            e.grid(row=0, column=2 * i + 1, padx=(2, 8))
            e.bind("<Return>", lambda _e: self.request_solve())
        ttk.Label(left, text="Exclude the liquid–vapour zone;\nthe model does "
                            "not represent it.", foreground="#52514e",
                  font=("", 8)).grid(row=4, column=0, columnspan=2, sticky="w",
                                     pady=(0, 8))

        row = 5
        for name in SLIDER_ORDER:
            spec = PARAMETERS[name]
            ttk.Label(left, text=spec["label"]).grid(row=row, column=0,
                                                     columnspan=2, sticky="w")
            lo, hi = spec["bounds"]
            if spec["log"]:
                lo, hi = np.log10(lo), np.log10(hi)
            s = ttk.Scale(left, from_=lo, to=hi, variable=self.values[name],
                          orient="horizontal", length=240,
                          command=lambda _v: self._sync_readouts())
            s.grid(row=row + 1, column=0, sticky="ew")
            s.bind("<ButtonRelease-1>", lambda _e: self.request_solve())
            ttk.Label(left, textvariable=self.readout[name], width=11,
                      anchor="e").grid(row=row + 1, column=1, sticky="e")
            row += 2

        ttk.Separator(left, orient="horizontal").grid(
            row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1
        for text, cmd in (("Auto-fit A_NH3 + n_H2", lambda: self.auto_fit(("A_NH3", "n_H2"))),
                          ("Auto-fit A_NH3 only", lambda: self.auto_fit(("A_NH3",))),
                          ("Reset to Kesten defaults", self.reset)):
            ttk.Button(left, text=text, command=cmd).grid(
                row=row, column=0, columnspan=2, sticky="ew", pady=1)
            row += 1

        ttk.Label(left, textvariable=self.metrics, justify="left",
                  font=("Consolas", 9)).grid(row=row, column=0, columnspan=2,
                                             sticky="w", pady=(10, 0))
        row += 1
        ttk.Label(left, textvariable=self.status, foreground="#2a78d6").grid(
            row=row, column=0, columnspan=2, sticky="w", pady=(6, 0))

        fig = Figure(figsize=(7.6, 5.6), facecolor=plots.SURFACE)
        self.ax_T = fig.add_subplot(211)
        self.ax_X = fig.add_subplot(212, sharex=self.ax_T)
        fig.tight_layout()
        self.fig = fig
        self.canvas = FigureCanvasTkAgg(fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def _sync_readouts(self):
        for n in SLIDER_ORDER:
            v = self._from_slider(n, self.values[n].get())
            self.readout[n].set(f"{v:.4g}" if not PARAMETERS[n]["log"]
                                else f"{v:.3e}")

    # -- solving on a worker thread ---------------------------------------
    def request_solve(self):
        """Queue a solve. tkinter is single threaded, so this must not block."""
        self.case.z_window = (self.zmin.get() * 1e-3, self.zmax.get() * 1e-3)
        if self.busy:
            self.pending = True
            return
        self.busy = True
        self.status.set("solving…")
        params = self.params()
        threading.Thread(target=self._worker, args=(params,), daemon=True).start()

    def _worker(self, params):
        try:
            errs, sol = self.case.errors(params)
            self.queue.put(("solved", params, errs, sol))
        except Exception as exc:                       # noqa: BLE001
            self.queue.put(("error", str(exc)))

    def auto_fit(self, names):
        if self.busy:
            return
        self.busy = True
        self.status.set(f"fitting {', '.join(names)} …")

        def run():
            try:
                def progress(n, _p, rms):
                    self.queue.put(("progress", n, rms))
                res = fit(self.case, fit_names=names, start=self.params(),
                          callback=progress)
                self.queue.put(("fitted", res))
            except Exception as exc:                   # noqa: BLE001
                self.queue.put(("error", str(exc)))

        threading.Thread(target=run, daemon=True).start()

    def _poll(self):
        try:
            while True:
                msg = self.queue.get_nowait()
                kind = msg[0]
                if kind == "solved":
                    _, params, errs, sol = msg
                    self.solution = sol
                    self._redraw(sol, errs, params)
                    self.status.set("ready")
                    self.busy = False
                    if self.pending:
                        self.pending = False
                        self.request_solve()
                elif kind == "progress":
                    self.status.set(f"fitting… evaluation {msg[1]}, "
                                    f"rms {msg[2]:.4f}")
                elif kind == "fitted":
                    res = msg[1]
                    for n in SLIDER_ORDER:
                        if n in res["params"]:
                            self.values[n].set(self._to_slider(n, res["params"][n]))
                    self._sync_readouts()
                    self.solution = res["solution"]
                    self._redraw(res["solution"], res["errors"], res["params"],
                                 sensitivity=res["sensitivity"])
                    self.status.set(f"fit complete ({res['nfev']} evaluations)")
                    self.busy = False
                elif kind == "error":
                    self.status.set(f"error: {msg[1]}")
                    self.busy = False
        except queue.Empty:
            pass
        self.root.after(60, self._poll)

    # -- drawing ----------------------------------------------------------
    def _redraw(self, sol, errs, params, sensitivity=None):
        d = self.case.data
        m = self.case.mask()
        for ax in (self.ax_T, self.ax_X):
            ax.clear()
            ax.grid(True, color=plots.GRID, linewidth=0.7)
            ax.set_axisbelow(True)

        self.ax_T.plot(sol.z * 1e3, sol.T_gas, color=plots.SERIES[0], lw=2,
                       label="Model, gas")
        self.ax_T.plot(sol.z * 1e3, sol.T_solid, color=plots.SERIES[1], lw=1.2,
                       ls=(0, (4, 3)), label="Model, catalyst")
        if d.T is not None:
            self.ax_T.plot(d.z[m] * 1e3, d.T[m], "o", ms=6,
                           color=plots.INK, label="Measured (fitted)")
            if (~m).any():
                self.ax_T.plot(d.z[~m] * 1e3, d.T[~m], "x", ms=6,
                               color=plots.MUTED, label="Excluded")
        self.ax_T.set_ylabel("Temperature [K]")
        self.ax_T.legend(loc="lower right", fontsize=8, frameon=False)

        self.ax_X.plot(sol.z * 1e3, sol.dissociation, color=plots.SERIES[2], lw=2,
                       label="Model")
        if d.X is not None:
            self.ax_X.plot(d.z[m] * 1e3, d.X[m], "o", ms=6, color=plots.INK,
                           label="Measured (fitted)")
            if (~m).any():
                self.ax_X.plot(d.z[~m] * 1e3, d.X[~m], "x", ms=6,
                               color=plots.MUTED, label="Excluded")
        self.ax_X.set_ylabel("Dissociation, X")
        self.ax_X.set_xlabel("Axial position, z [mm]")
        self.ax_X.set_ylim(0, 1)
        self.ax_X.legend(loc="lower right", fontsize=8, frameon=False)

        for ax in (self.ax_T, self.ax_X):
            ax.axvspan(ax.get_xlim()[0], self.zmin.get(), color=plots.GRID, alpha=0.5)

        # Only report channels that are actually weighted into the objective.
        # Pressure is carried for plotting but not fitted, and showing its
        # residual next to the others invites reading it as a fit quality.
        units = {"T": "K", "X": "", "p": "Pa"}
        lines = [f"{'RMS ' + k:<9}{v:10.4g} {units.get(k, '')}"
                 for k, v in errs.items() if k in ("T", "X")]
        lines.append(f"{'points':<9}{int(m.sum()):10d} of {len(d)}")
        if int(m.sum()) < 4:
            lines.append("few points: expect a")
            lines.append("poorly determined fit")
        if sensitivity:
            lines.append("")
            lines.append("sensitivity (10% change):")
            for k, v in sensitivity.items():
                flag = "  <- unconstrained" if v < 0.02 else ""
                lines.append(f"  {k:<8}{v:8.4f}{flag}")
        self.metrics.set("\n".join(lines))
        self.fig.tight_layout()
        self.canvas.draw_idle()

    # -- actions ----------------------------------------------------------
    def load_csv(self):
        path = filedialog.askopenfilename(
            title="Measured bed data",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.case.data = BedData.from_csv(path)
        except Exception as exc:                        # noqa: BLE001
            self.status.set(f"could not load: {exc}")
            return
        self.zmax.set(self.case.data.z[-1] * 1e3)
        self.root.title(f"Catalyst bed kinetics — {self.case.data.name}")
        self.request_solve()

    def reset(self):
        for n in SLIDER_ORDER:
            self.values[n].set(self._to_slider(n, PARAMETERS[n]["default"]))
        self._sync_readouts()
        self.request_solve()


def build_case(args):
    if args.data is None:
        return FitCase.kesten_demo()
    data = BedData.from_csv(args.data)
    bed = CatalystBed.uniform(
        args.bed_dia * 1e-3, args.bed_length * 1e-3,
        PackedSpheres(mesh=(args.mesh_coarse, args.mesh_fine), eps=args.voidage))
    return FitCase(bed, G=bed.mass_flux(args.mdot * 1e-3),
                   p_inlet=args.p_inlet * 1e5, data=data, T_feed=args.t_feed)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data', default=None, help='measured CSV; omit for the Kesten demo')
    ap.add_argument('--bed-dia', type=float, default=8.3, metavar='MM')
    ap.add_argument('--bed-length', type=float, default=20.0, metavar='MM')
    ap.add_argument('--mdot', type=float, default=0.0443, metavar='G_S')
    ap.add_argument('--p-inlet', type=float, default=8.46, metavar='BAR')
    ap.add_argument('--t-feed', type=float, default=298.15, metavar='K')
    ap.add_argument('--mesh-coarse', type=int, default=25)
    ap.add_argument('--mesh-fine', type=int, default=30)
    ap.add_argument('--voidage', type=float, default=0.38)
    args = ap.parse_args()

    np.seterr(all='ignore')
    root = tk.Tk()
    FitApp(root, build_case(args))
    root.mainloop()


if __name__ == '__main__':
    main()
