"""Plot held-out rates from the archived benchmark summaries."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

DATA = Path(__file__).resolve().parent / "data"
METHODS = ("Crosslabel", "Marked pairs")
COLORS = ("#087F82", "#645190")
MARKERS = ("o", "s")
LABELS = (
    "Original matched-pair setting",
    "Lower density, balanced codes",
    "Low retention, skewed codes",
    "Varying intensity, skewed",
)


def load(path):
    return json.loads(path.read_text())


def original_row(summary, name):
    source = summary[name]
    return (
        (source["first_alerts"], 100, source["first_ci"]),
        (source["second_alerts"], 100, source["second_ci"]),
    )


def new_row(source):
    return tuple(
        (source[f"{method}_alerts"], source["n"], source[f"{method}_ci"])
        for method in ("nested", "pair")
    )


def draw_panel(ax, rows, title, right, ticks, show_labels):
    ax.set_title(title, loc="left", fontsize=9.3, fontweight="bold", pad=10)
    ax.set_ylim(-.46, 3.46)
    ax.set_xlim(-.024 if right == 1 else -.009, right)
    ax.set_yticks(range(4), LABELS if show_labels else [""] * 4)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0, labelsize=8.2)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", color="#DFE6EC", linewidth=.7)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#A9B5BF")
    ax.axhline(.5, color="#C8D2DC", linewidth=.8)

    for row, pair in enumerate(rows):
        if pair is None:
            ax.text(.01, row, "Not evaluated", color="#526575",
                    fontsize=7.2, va="center")
            continue
        for (count, n, (lo, hi)), color, marker, offset in zip(
            pair, COLORS, MARKERS, (-.135, .135)
        ):
            rate = count / n
            ax.errorbar(rate, row + offset, xerr=[[max(0., rate-lo)], [max(0., hi-rate)]],
                        marker=marker, color=color, markersize=4, linewidth=1, capsize=2)
            ax.text(1.03, row + offset, f"{count}/{n}", transform=ax.get_yaxis_transform(),
                    color=color, fontsize=7.7, va="center", clip_on=False)


def main(destination):
    previous = load(DATA / "marked_pair/SUMMARY.json")
    recent = load(DATA / "robustness_reference/SUMMARY.json")["settings"]
    settings = [
        recent["balanced_lower_density"],
        recent["low_retention_skewed_codes"],
        recent["varying_intensity_skewed_codes"],
    ]
    detections = [original_row(previous, "matched_alternative")] + [
        new_row(setting["None"]) for setting in settings
    ]
    on_grid = [original_row(previous, "matched_adequate")] + [
        new_row(setting["on_grid"]) for setting in settings
    ]
    off_grid = [None] + [new_row(setting["off_grid"]) for setting in settings]

    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": "#203344",
                         "axes.labelcolor": "#203344", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(7.08, 2.44), facecolor="white")
    ax1 = fig.add_axes([.258, .26, .207, .57])
    ax2 = fig.add_axes([.55, .26, .156, .57])
    ax3 = fig.add_axes([.788, .26, .154, .57])
    draw_panel(ax1, detections, "A  Trimer detections", 1., [0, .5, 1], True)
    draw_panel(ax2, on_grid, "B  On-grid alerts", .24, [0, .1, .2], False)
    draw_panel(ax3, off_grid, "C  Off-grid alerts", .24, [0, .1, .2], False)
    handles = [
        Line2D([], [], color=color, marker=marker, linestyle="-", markersize=4,
               linewidth=1, label=method)
        for method, color, marker in zip(METHODS, COLORS, MARKERS)
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(.575, .015), fontsize=8)
    fig.savefig(destination, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)

    rows = [
        json.loads(line)
        for line in (DATA / "marked_pair/NPC_TEMPLATES.jsonl")
        .read_text().splitlines()
    ]
    controls = {row["seed"]: row["ci"][0] for row in rows if row["adequate"]}
    alternatives = {row["seed"]: row["ci"][0] for row in rows if not row["adequate"]}
    assert len(rows) == 200 and len(controls) == len(alternatives) == 100
    assert controls.keys() == alternatives.keys()
    fig, ax = plt.subplots(figsize=(3.4, 1.87), layout="constrained")
    ax.axvspan(0, .12, color="#F0E9F5")
    ax.axvline(0, color="#645190", linewidth=.8)
    ax.scatter(list(controls.values()), [alternatives[seed] for seed in controls],
               s=13, color="#087F82", alpha=.78, edgecolors="none")
    ax.set(xlim=(-.9, .12), ylim=(3.2, 6.1), yticks=[3.5, 4.5, 5.5],
           xlabel="Lower endpoint, two-code control",
           ylabel="Lower endpoint, alternative")
    ax.tick_params(labelsize=8)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="y", color="#DFE6EC", linewidth=.7)
    ax.set_axisbelow(True)
    fig.savefig(destination.with_name("viability_templates.pdf"),
                metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    main(parser.parse_args().destination)
