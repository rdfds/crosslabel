"""Render the three-code field-check schematic."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch


BLUE, ORANGE, PURPLE, INK, FAINT = "#076c83", "#b6541e", "#7555a2", "#172a38", "#dde5ea"
COLORS = [BLUE, ORANGE, PURPLE]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK, "axes.edgecolor": "#94a4ae", "axes.titleweight": "bold",
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42,
})



def panel(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", fontsize=10, pad=9)


def workflow(output):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.28), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1, 1.34]})
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    for ax, codes, label, letter, title in [
        (axes[0], [0, 1, 1], "Two distinct codes", "A", "Restricted component"),
        (axes[1], [0, 1, 2], "Three distinct codes", "B", "Three-code example"),
    ]:
        panel(ax, letter, title)
        points = [(0.47, .48), (.25, .68), (.70, .63)]
        ax.add_patch(Circle(points[0], .31, ls="--", lw=1.15, fill=False, color="#8a9ca8"))
        for x, y in points[1:]:
            ax.plot([points[0][0], x], [points[0][1], y], color="#8a9ca8", lw=1.1)
        for (x, y), code in zip(points, codes):
            ax.scatter(x, y, s=110, color=COLORS[code], edgecolor="white", linewidth=.8, zorder=3)
        ax.text(.47, .06, label, ha="center", weight="bold", fontsize=8.7)
        ax.text(.47, .235, "anchor radius", ha="center", color="#526772", fontsize=8)
    ax = axes[2]
    panel(ax, "C", "Ungrouped field")
    for y, heading, detail, color in [
        (.76, "Estimate connected signal", "Local counts and background correction", INK),
        (.47, "Check support and uncertainty", "Spatial scale and calibrated limits", INK),
        (.18, "If lower limit exceeds zero", "Review assumptions or expand library", BLUE),
    ]:
        ax.add_patch(FancyBboxPatch((.015, y-.105), .96, .21, boxstyle="round,pad=.012",
                                    facecolor="#f3f7f9", edgecolor=FAINT, lw=.8))
        ax.text(.06, y+.038, heading, fontsize=8.6, weight="bold", color=color, va="center")
        ax.text(.06, y-.044, detail, fontsize=7.3, va="center")
    for y in (.618, .328):
        ax.annotate("", xy=(.50, y-.041), xytext=(.50, y+.041),
                    arrowprops=dict(arrowstyle="-|>", color="#899daa", lw=1))
    fig.savefig(output / "fig1_workflow.pdf", facecolor="white",
                metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    output = parser.parse_args().output
    output.mkdir()
    workflow(output)
