"""Report figures from the SWMM results (print, 300 dpi). Usage: python make_figures.py"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "figures"
FIG.mkdir(exist_ok=True)

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
F2R, SG, NONE = "#2a78d6", "#eb6834", "#898781"   # validated pair; "no unit" is a neutral reference

plt.rcParams.update({
    "font.family": "Calibri", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS, "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False, "savefig.dpi": 300,
})


def hydrograph():
    h = pd.read_csv(RES / "events" / "event_hydrographs.csv", index_col=0, parse_dates=True)
    t0 = h.index[0]
    hrs = (h.index - t0).total_seconds() / 3600
    fig, ax = plt.subplots(figsize=(6.3, 2.9))
    series = [("mar2026_mbz|none", "No unit", NONE, (0, (4, 2))),
              ("mar2026_mbz|singapore", "Singapore-style rain garden", SG, "-"),
              ("mar2026_mbz|f2r", "Flood2Resource", F2R, "-")]
    for col, label, c, ls in series:
        ax.plot(hrs, h[col], color=c, lw=2, ls=ls, label=label, solid_capstyle="round")
    ax.set_xlim(8, 20)
    ax.set_ylim(0, 60)
    ax.set_xlabel("Time from start of storm (h)")
    ax.set_ylabel("Flow to municipal drain (L/s)")
    peak_none, peak_f2r = h["mar2026_mbz|none"].max(), h["mar2026_mbz|f2r"].max()
    ax.annotate(f"{peak_none:.0f} L/s", xy=(12.05, peak_none), xytext=(13.0, 55), color=INK2,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.annotate(f"{peak_f2r:.0f} L/s once the tank is full\n(peak cut by 81%)", xy=(13.9, peak_f2r),
                xytext=(15.0, 24), color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=8.5, handlelength=2.6)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_hydrograph_mar2026_mbz.png")
    plt.close(fig)


def peak_bars():
    r = pd.read_csv(RES / "events" / "event_results.csv")
    # two groups: ordinary storms (what the Singapore design is built for) and desert flash floods
    order = [("typical_10", "10 mm storm"), ("typical_20", "20 mm storm"), ("typical_30", "30 mm storm"),
             None,
             ("burst_30min", "30 mm in 30 minutes"), ("mar2026_mbz", "March 2026, MBZ City"),
             ("mar2026_wathba", "March 2026, Al Wathba")]
    fig, ax = plt.subplots(figsize=(6.3, 3.3))
    h = 0.36
    rows = [i for i, o in enumerate(order) if o]
    for k, (cfg, label, c) in enumerate([("singapore", "Original Singapore design", SG), ("f2r", "Flood2Resource", F2R)]):
        vals = [r[(r.event == order[i][0]) & (r.config == cfg)].peak_reduction_pct.iloc[0] for i in rows]
        pos = [i + (k - 0.5) * (h + 0.04) for i in rows]
        ax.barh(pos, vals, height=h, color=c, label=label)
        for p, v in zip(pos, vals):
            ax.text(max(v, 0) + 1.5, p, f"{v:.0f}%", va="center", fontsize=8, color=INK2)
    ax.set_yticks(rows, [order[i][1] for i in rows], color=INK2)
    ax.invert_yaxis()
    ax.axhline(3, color=AXIS, lw=0.8)
    ax.text(-0.02, -0.75, "Ordinary storms", transform=ax.get_yaxis_transform(), ha="right", fontsize=8.5,
            fontweight="bold", color=INK)
    ax.text(-0.02, 3.25, "Desert flash floods", transform=ax.get_yaxis_transform(), ha="right", fontsize=8.5,
            fontweight="bold", color=INK)
    ax.set_ylim(len(order) - 0.4, -1.1)
    ax.set_xlim(0, 110)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Reduction in peak flow to the drain (%)")
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_peak_reduction_by_storm.png")
    plt.close(fig)


def combined():
    """(a) hydrograph and (b) peak reduction by storm side by side, one shared legend, sized for 16 cm width."""
    h = pd.read_csv(RES / "events" / "event_hydrographs.csv", index_col=0, parse_dates=True)
    r = pd.read_csv(RES / "events" / "event_results.csv")
    hrs = (h.index - h.index[0]).total_seconds() / 3600
    with plt.rc_context({"font.size": 7.5}):
        fig, (a, b) = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"width_ratios": [1, 1.12]})
        # (a) hydrograph
        for col, label, c, ls in [("mar2026_mbz|none", "No unit", NONE, (0, (4, 2))),
                                  ("mar2026_mbz|singapore", "Original Singapore design", SG, "-"),
                                  ("mar2026_mbz|f2r", "Flood2Resource", F2R, "-")]:
            a.plot(hrs, h[col], color=c, lw=1.6, ls=ls, label=label, solid_capstyle="round")
        a.set_xlim(8, 20); a.set_ylim(0, 58)
        a.set_xlabel("Hours from start of storm"); a.set_ylabel("Flow to municipal drain (L/s)")
        pf = h["mar2026_mbz|f2r"].max()
        a.annotate(f"{pf:.0f} L/s after the tank fills\n(peak cut by 81%)", xy=(13.9, pf), xytext=(14.6, 24),
                   color=INK2, fontsize=7, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
        a.set_title("(a) March 2026 storm, MBZ City street", fontsize=8, loc="left", color=INK)
        # (b) peak reduction
        order = [("typical_10", "10 mm"), ("typical_20", "20 mm"), ("typical_30", "30 mm"), None,
                 ("burst_30min", "30 mm in 30 min"), ("mar2026_mbz", "Mar 2026, MBZ City"),
                 ("mar2026_wathba", "Mar 2026, Al Wathba")]
        rows = [i for i, o in enumerate(order) if o]
        hh = 0.38
        for k, (cfg, c) in enumerate([("singapore", SG), ("f2r", F2R)]):
            vals = [r[(r.event == order[i][0]) & (r.config == cfg)].peak_reduction_pct.iloc[0] for i in rows]
            pos = [i + (k - 0.5) * (hh + 0.03) for i in rows]
            b.barh(pos, vals, height=hh, color=c)
            for p, v in zip(pos, vals):
                b.text(max(v, 0) + 1.5, p, f"{v:.0f}%", va="center", fontsize=6.5, color=INK2)
        b.set_yticks(rows, [order[i][1] for i in rows], color=INK2)
        b.invert_yaxis(); b.set_ylim(len(order) - 0.45, -1.2)
        b.axhline(3, color=AXIS, lw=0.7)
        b.text(108, -0.85, "Ordinary storms", ha="right", fontsize=7, fontweight="bold", color=INK)
        b.text(108, 3.35, "Desert flash floods", ha="right", fontsize=7, fontweight="bold", color=INK)
        b.set_xlim(0, 112); b.set_xticks([0, 25, 50, 75, 100])
        b.grid(axis="y", visible=False); b.tick_params(axis="y", length=0)
        b.set_xlabel("Reduction in peak flow to the drain (%)")
        b.set_title("(b) Peak reduction by storm", fontsize=8, loc="left", color=INK)
        handles, labels = a.get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=7.5)
        fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=1.5)
        fig.savefig(FIG / "fig_results_combined.png")
        plt.close(fig)


if __name__ == "__main__":
    hydrograph()
    peak_bars()
    combined()
    print("figures written to", FIG)
