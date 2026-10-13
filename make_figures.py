"""Figures for the README from results/docking_summary.csv and results/docking_scores_all_seeds.csv."""
import matplotlib.pyplot as plt
import pandas as pd

COLORS = {"control": "#2a78d6", "candidate": "#eb6834", "background": "#a3a29c"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.family": "Helvetica", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK})

summary = pd.read_csv("results/docking_summary.csv")
seeds = pd.read_csv("results/docking_scores_all_seeds.csv")
bg = summary[summary.group == "background"]

# --- Figure 1: scores vs the random-drug background ---
ranked = summary[summary.group != "background"].sort_values("mean", ascending=False)
fig, ax = plt.subplots(figsize=(8, 5.6))
lo, med, hi = bg["mean"].quantile(0.1), bg["mean"].median(), bg["mean"].quantile(0.9)
ax.axvspan(lo, hi, color=GRID, zorder=0)
ax.axvline(med, color=COLORS["background"], lw=1.5, ls="--", zorder=1)
ax.text(med, len(ranked) + 0.35, f"random approved drugs (n={len(bg)})\nmedian {med:.1f}, 10th-90th pct shaded",
        ha="center", va="bottom", fontsize=8.5, color=MUTED)
for y, (_, r) in enumerate(ranked.iterrows()):
    pts = seeds[seeds.compound == r.compound].score
    ax.plot([pts.min(), pts.max()], [y, y], color=COLORS[r.group], lw=2, alpha=0.45, zorder=2)
    ax.scatter(r["mean"], y, s=64, color=COLORS[r.group], edgecolor="white", linewidth=2, zorder=3)
    ax.text(r["mean"] - 0.18, y, f"{r['mean']:.1f}", va="center", ha="left", fontsize=8.5, color=INK)
ax.set_yticks(range(len(ranked)), ranked.compound)
ax.set_xlabel("AutoDock Vina score, kcal/mol (more negative = stronger predicted binding)")
ax.invert_xaxis()
ax.set_ylim(-0.7, len(ranked) + 1.3)
ax.grid(axis="x", color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
handles = [plt.Line2D([], [], marker="o", ls="", ms=8, color=COLORS[g], label=l)
           for g, l in [("control", "Approved IDH1 inhibitors (positive controls)"),
                        ("candidate", "Screening candidates")]]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
ax.set_title("Docking into the IDH1 R132H allosteric pocket (PDB 8T7O)\n"
             "Dot = mean of 3 seeds, line = range", loc="left", fontsize=11, color=INK)
fig.tight_layout()
fig.savefig("results/docking_scores.png", dpi=200, facecolor="white")

# --- Figure 2: size bias ---
fig, ax = plt.subplots(figsize=(7, 4.6))
for g, label in [("background", "Random approved drugs"), ("candidate", "Screening candidates"),
                 ("control", "Positive controls")]:
    d = summary[summary.group == g]
    ax.scatter(d.heavy_atoms, d["mean"], s=56, color=COLORS[g], edgecolor="white", linewidth=2,
               label=label, zorder=3)
for _, r in summary[summary.group != "background"].iterrows():
    if r.group == "control" or r["mean"] <= ranked["mean"].min() + 2:
        ax.annotate(r.compound, (r.heavy_atoms, r["mean"]), xytext=(6, 2),
                    textcoords="offset points", fontsize=8.5, color=INK)
r = summary.heavy_atoms.corr(-summary["mean"])  # size vs predicted affinity
ax.set_xlabel("Heavy atoms (molecule size)")
ax.set_ylabel("Vina score, kcal/mol")
ax.invert_yaxis()
ax.grid(color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.set_title(f"Bigger molecules score better\nPearson r = {r:.2f} between size and predicted affinity, "
             f"n = {len(summary)} drugs", loc="left", fontsize=11, color=INK)
fig.tight_layout()
fig.savefig("results/size_bias.png", dpi=200, facecolor="white")
print(f"Saved figures. Size-score correlation r = {r:.2f}")
