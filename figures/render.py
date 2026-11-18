"""Build the Crosslabel manuscript artwork at its final print size."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator, PercentFormatter

import style as s

METHODS = {
    "crosslabel_lower": (s.CROSS, "-"),
    "crosslabel_point": (s.CROSS, "--"),
    "dtc_connected": (s.TRIPLE, "-."),
    "marked_pair": (s.PAIR, "-"),
}


def circle(ax, center, radius, linestyle="-", color=s.REF):
    ax.add_patch(Circle(center, radius, fill=False, edgecolor=color,
                        linewidth=.6, linestyle=linestyle, zorder=4))


def record_map(ax, xy, codes, extent, size, *, alpha=.88):
    ax.set(xlim=(-extent, extent), ylim=(-extent, extent))
    ax.set_axis_off()
    xy, codes = np.asarray(xy), np.asarray(codes)
    for code, color in enumerate(s.ROUNDS):
        subset = xy[codes == code]
        ax.scatter(*subset.T, color=color, s=size, alpha=alpha, linewidths=0)
    return ax


def mark_witness(ax, points, *, size):
    anchor, *neighbors = np.asarray(points)
    ax.scatter(*np.asarray(neighbors).T, s=size, facecolors="none", edgecolors=s.INK,
               linewidths=.65, zorder=6)
    ax.scatter(*anchor, s=size * 1.25, facecolors="none", edgecolors=s.INK,
               linewidths=.8, marker="s", zorder=6)


def code_key(fig, box, x, y, *, separation=8, vertical=False):
    from matplotlib.lines import Line2D
    for code, color in enumerate(s.ROUNDS):
        start = x if vertical else x + code * separation
        height = y + code * separation if vertical else y
        fig.add_artist(Line2D(*([v] for v in s.mm(fig, box, start, height)),
                              transform=fig.transFigure, marker="o",
                              color=color, markersize=2.8, linestyle="none"))
        s.panel_text(fig, box, start + 1.45, height, f"R{code + 1}", va="center")


def arrow(fig, box, x0, y0, x1, y1):
    fig.add_artist(FancyArrowPatch(s.mm(fig, box, x0, y0), s.mm(fig, box, x1, y1),
                                  transform=fig.transFigure, arrowstyle="->",
                                  mutation_scale=6, lw=.6, color=s.REF,
                                  shrinkA=0, shrinkB=0))


def workflow(output, model, compact):
    if compact:
        fig, (left, right) = s.canvas(178, 37, [89, 82], [3])
        component_size, positions, y, map_x, map_size = 27.5, (16, 51), 4.5, 19, 28
    else:
        fig, (left, right, rules) = s.canvas(178, 55, [58, 40, 70], [3, 3])
        component_size, positions, y, map_x, map_size = 26, (1, 29), 9, 2, 36
    s.panel_letter(fig, left, "a")
    s.panel_letter(fig, right, "b")
    limit = max(39, 2 + max(abs(v) for item in model["components"] for row in item["xy"] for v in row))
    for i, (label, view) in enumerate(zip(("Dimer", "Trimer"), model["components"])):
        ax = s.axis(fig, left, positions[i], y, component_size, component_size, aspect="equal")
        record_map(ax, view["xy"], view["codes"], limit, 4.8)
        circle(ax, (0, 0), model["radius_nm"])
        ax.scatter(0, 0, s=14, marker="s", facecolors="none",
                   edgecolors=s.INK, linewidths=.7, zorder=6)
        s.panel_text(fig, left, positions[i] + component_size / 2,
                     y - 1.9, label, ha="center", va="bottom")
        if i == 0:
            scale = model["scale_bar_nm"]
            s.bar(ax, limit - scale - 2, -limit + 2, scale, f"{scale} nm")

    crop = model["crop"]
    s.panel_text(fig, right, map_x + map_size / 2, 2.5 if compact else 6,
                 "Monomer–dimer field", ha="center", va="bottom")
    ax = s.axis(fig, right, map_x, 4 if compact else 8, map_size, map_size, aspect="equal")
    record_map(ax, crop["xy"], crop["codes"], model["crop_nm"] / 2, 2.2, alpha=.68)
    for radius in model["annulus_nm"]:
        circle(ax, (0, 0), radius, (0, (2, 2)))
    circle(ax, (0, 0), model["radius_nm"])
    mark_witness(ax, crop["witness"], size=7)
    edge = model["crop_nm"] / 2
    for x in (-edge, edge):
        for y_corner in (-edge, edge):
            ax.plot([x, x * .93], [y_corner, y_corner], color=s.REF, lw=.6, clip_on=False)
            ax.plot([x, x], [y_corner, y_corner * .93], color=s.REF, lw=.6, clip_on=False)
    if compact:
        code_key(fig, right, 59, 10, separation=5.8, vertical=True)
    else:
        code_key(fig, left, 10, 43, separation=10)
        s.panel_letter(fig, rules, "c")
        for x, label in ((6.5, "Localizations\nwith round\ncodes"),
                         (24.2, "Connected\nthree-code\ncontrast"),
                         (39.3, "Spatial\ninterval")):
            s.panel_text(fig, rules, x, 16, label, ha="center", va="center")
        arrow(fig, rules, 14, 16, 17.2, 16)
        arrow(fig, rules, 31.3, 16, 33.8, 16)
        arrow(fig, rules, 44.1, 15, 47, 10)
        arrow(fig, rules, 44.1, 17, 47, 22)
        s.panel_text(fig, rules, 48, 9, "Interval rule\nLower limit > 0", va="center")
        s.panel_text(fig, rules, 48, 23, "Reference rule\nRank against\nsimulated adequate\nfields", va="center")
        # A single junction makes both outcomes available to either rule.
        arrow(fig, rules, 57, 31, 57, 34)
        s.panel_text(fig, rules, 56, 35, "Flag: revisit library\nNo flag: retain library\nfor fitting",
                     ha="center", va="top")
    s.save(fig, output, "crosslabel_fig1")


def matched(output, panels):
    fig, boxes = s.canvas(86, 62, [38.5, 38.5], [5])
    rng = np.random.default_rng(panels[0]["groups"][0]["fields"][0]["seed"])
    jitter = [rng.uniform(-.2, .2, group["n"]) for group in panels[0]["groups"]]
    for i, (box, panel) in enumerate(zip(boxes, panels)):
        s.panel_letter(fig, box, "ab"[i])
        ax = s.axis(fig, box, 12.6, 12, 24.6, 32.7)
        scores = np.asarray([r["score"] for group in panel["groups"] for r in group["fields"]])
        low, high = min(scores.min(), panel["cutoff"]), max(scores.max(), panel["cutoff"])
        pad = (high - low) * .075
        ax.set(xlim=(-.55, 1.55), ylim=(low - pad, high + pad),
               xticks=[0, 1], xticklabels=["Adequate", "Trimer"])
        if i == 0:
            ax.set_ylabel("Marked-pair squared residual", labelpad=2.0)
            ax.set_yticks([0, 100, 200, 300] if high + pad > 300 else [50, 100, 150, 200])
            ax.set_ylim(low - pad, high + pad)
        else:
            ax.set_ylabel("Lower endpoint\n(pairs per anchor)", labelpad=2.0)
            ticks = MaxNLocator(nbins=4).tick_values(low - pad, high + pad)
            ax.set_yticks([value for value in ticks if low - pad < value < high + pad])
            ax.yaxis.set_major_formatter(FuncFormatter(
                lambda value, _: f"{value:.2f}".replace("-", "−")
            ))
        for j, (group, color) in enumerate(zip(panel["groups"], (s.ADEQUATE, s.ALTERNATIVE))):
            ax.scatter(j + jitter[j], [field["score"] for field in group["fields"]],
                       color=color, alpha=.6, s=6.25, linewidths=0)
            ax.text(j, 1.072, f'{group["flags"]}/{group["n"]}\nflagged',
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom")
        ax.axhline(panel["cutoff"], color=s.REF, lw=.6, linestyle="--", zorder=0)
    s.save(fig, output, "crosslabel_fig2")


def coverage_panel(fig, box, model, *, standalone=False):
    rows = model["rows"]
    y = np.r_[np.arange(6), np.arange(7, 13)]
    if standalone:
        ax = s.axis(fig, box, 42, 5, 73, 37)
    else:
        ax = s.axis(fig, box, 27.5, 5, 24, 40)
    ax.set(xlim=(.84, 1.003), ylim=(12.62, -.61),
           xticks=[.85, .90, .95, 1.0],
           yticks=y, yticklabels=[r["label"] for r in rows])
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_xlabel("Interval coverage (%)", labelpad=2.3)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["left"].set_visible(False)
    lo, hi, n = model["acceptance_counts"]
    ax.add_patch(Rectangle((lo / n, -.46), (hi - lo) / n, 5.92,
                           facecolor="#eeeeee", edgecolor="none", zorder=0))
    ax.axvline(model["nominal"], color=s.REF, lw=.6, linestyle="--", zorder=1)
    for yy, row in zip(y, rows):
        value = row["rate"]
        left, right = row["ci"]
        ax.errorbar(value, yy, xerr=[[value - left], [right - value]],
                    color=s.CROSS, fmt="o", mfc=s.CROSS if row["primary"] else "white",
                    markersize=2.8, lw=1, capsize=1.25, capthick=.65, zorder=4)


def roc(output, panels, coverage, separate):
    if separate:
        fig, boxes = s.canvas(178, 46, [58, 53, 53], [5, 5])
        ax_w, ax_h = 45, 22.5
    else:
        fig, boxes = s.canvas(178, 58, [41, 34, 34, 56], [3, 3, 3])
        ax_w, ax_h = 27, 34
    legend_methods = (
        ("crosslabel_lower", "Crosslabel lower"),
        ("crosslabel_point", "Crosslabel point"),
        ("dtc_connected", "Connected triple"),
        ("marked_pair", "Marked pair"),
    )
    for index, (box, panel) in enumerate(zip(boxes, panels)):
        s.panel_letter(fig, box, "abc"[index])
        ax_x = 9 if index == 0 else (4 if separate else 3)
        ax = s.axis(fig, box, ax_x, 7, ax_w, ax_h)
        ax.set(xlim=(-.65, 20.55), ylim=(-2, 82.1),
               xticks=[0, 5, 10, 15, 20] if separate else [0, 5, 10, 20],
               yticks=[0, 20, 40, 60, 80])
        ax.tick_params(axis="y", labelleft=index == 0)
        if index == 0:
            ax.set_ylabel("Detections (of 80)", labelpad=2)
        if index == 1:
            ax.set_xlabel("False alerts (of 100)", labelpad=2.5)
            legend_center = ax.get_position().x0 + ax.get_position().width / 2
        s.panel_text(fig, box, ax_x + ax_w / 2, 5.2, panel["label"], ha="center", va="bottom")
        ax.axvline(panel["comparison_alerts"], color=s.REF, lw=.6, linestyle=":", zorder=0)
        for method in ("marked_pair", "crosslabel_lower", "dtc_connected", "crosslabel_point"):
            color, linestyle = METHODS[method]
            points = np.array(panel["curves"][method])
            x, y = points.T
            if x[-1] < 20:
                x, y = np.r_[x, 20], np.r_[y, y[-1]]
            ax.plot(x, y, color=color, linestyle=linestyle, drawstyle="steps-post", zorder=2)
        # Place all nominal and local decisions at their recorded coordinates.
        for method in ("marked_pair", "crosslabel_lower", "dtc_connected", "crosslabel_point"):
            color, _ = METHODS[method]
            ax.plot(*panel["ranked"][method], color=color, marker="o", linestyle="none",
                    mfc="none", mew=.65, markersize=2.5 if method == "crosslabel_point" else 3.5,
                    zorder=5)
        for method in ("crosslabel_lower", "dtc_connected"):
            ax.plot(*panel["local"][method], color=METHODS[method][0], marker="s",
                    linestyle="none", mew=0, markersize=2.5, zorder=6)
    handles = [Line2D([], [], color=METHODS[method][0], linestyle=METHODS[method][1],
                      label=label) for method, label in legend_methods]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(legend_center, .025),
               ncols=4, frameon=False, handlelength=2.7, handletextpad=.5,
               columnspacing=1.1, borderaxespad=0, borderpad=0)
    if not separate:
        s.panel_letter(fig, boxes[3], "d")
        coverage_panel(fig, boxes[3], coverage)
    s.save(fig, output, "crosslabel_fig3")
    if separate:
        fig, (box,) = s.canvas(125, 56, [121], [])
        coverage_panel(fig, box, coverage, standalone=True)
        s.save(fig, output, "crosslabel_coverage")


def nup96(output, model):
    fig, (left, right) = s.canvas(178, 62, [73, 96], [5])
    s.panel_letter(fig, left, "a")
    s.panel_letter(fig, right, "b")
    crop_half = model["crop_nm"] / 2
    ax = s.axis(fig, left, 1.1, 6.5, 43, 43, aspect="equal")
    record_map(ax, model["crop_xy"], model["crop_codes"], crop_half, 4.2)
    inset_half = model["inset_nm"] / 2
    ax.add_patch(Rectangle((-inset_half, -inset_half), model["inset_nm"], model["inset_nm"],
                           fill=False, edgecolor=s.REF, lw=.6))
    scale = model["bars_nm"][0]
    s.bar(ax, -crop_half + 35, -crop_half + 45, scale, f"{scale} nm")
    # The locator shows the same coordinate window as the enlargement.
    inset = s.axis(fig, left, 48, 10, 23, 23, aspect="equal")
    record_map(inset, model["inset_xy"], model["inset_codes"], inset_half, 9)
    inset.add_patch(Rectangle((-inset_half, -inset_half), model["inset_nm"], model["inset_nm"],
                              fill=False, edgecolor=s.REF, lw=.6, clip_on=False))
    anchor = model["witness"][0]
    circle(inset, anchor, model["radius_nm"])
    mark_witness(inset, model["witness"], size=18)
    scale = model["bars_nm"][1]
    s.bar(inset, inset_half - scale - 17, -inset_half + 8, scale, f"{scale} nm")
    start = (inset_half / model["crop_nm"]) * 43
    ends = [s.mm(fig, left, 1.1 + 21.5 + start, 28.0),
            s.mm(fig, left, 45.2, 28.0), s.mm(fig, left, 48, 21.5)]
    fig.add_artist(Line2D(*zip(*ends), transform=fig.transFigure, color=s.REF, lw=.6, zorder=-1))
    code_key(fig, left, 48, 39.8, separation=6.5)

    ax = s.axis(fig, right, 24, 8, 67, 35)
    lower = np.asarray([field["lower"] for group in model["paired"] for field in group["fields"]])
    xmin, xmax = lower.min(), lower.max()
    pad = (xmax - xmin) * .06
    ax.set(xlim=(xmin - pad, xmax + pad), ylim=(-.62, 1.67),
           xticks=[-1, 0, 2, 4, 6], yticks=[1, 0],
           yticklabels=["Two-code control", "Alternative"])
    ax.tick_params(axis="y", length=0, pad=4)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Crosslabel lower endpoint (pairs per anchor)", labelpad=3)
    ax.axvline(model["zero_cutoff"], color=s.REF, lw=.6, zorder=0)
    rng = np.random.default_rng(model["paired"][0]["fields"][0]["seed"])
    jitter = rng.uniform(-.16, .16, model["paired"][0]["n"])
    for yy, group, color in zip((1, 0), model["paired"], (s.ADEQUATE, s.ALTERNATIVE)):
        values = np.array([field["lower"] for field in group["fields"]])
        label = f'{group["flags"]}/{group["n"]} flagged'
        ax.scatter(values, yy + jitter, color=color, s=6.25, alpha=.6, linewidths=0)
        if yy:
            positive = values > model["zero_cutoff"]
            ax.scatter(values[positive], (yy + jitter)[positive], color=s.INK,
                       s=4, linewidths=0, zorder=5)
            ax.annotate(label, (max(values), yy + max(jitter)),
                        xytext=(9, 5), textcoords="offset points", ha="left", va="bottom")
        else:
            ax.text(max(values), yy + .24, label, ha="right", va="bottom")
    s.save(fig, output, "crosslabel_fig4")


def rate_axes(fig, box, x, y, width, height, count, xmax, ticks):
    ax = s.axis(fig, box, x, y, width, height)
    ax.set(xlim=(-.006, xmax), ylim=(count - .5, -.5), yticks=[],
           xticks=ticks)
    ax.spines["left"].set_visible(False)
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    return ax


def rate_dot(ax, row, y, color, marker):
    value = row["count"] / row["n"]
    lo, hi = row["ci"]
    ax.errorbar(value, y, xerr=[[max(0., value - lo)], [max(0., hi - value)]],
                fmt=marker, color=color, mfc=color, mew=.7, lw=1,
                markersize=3, capsize=1.3, capthick=.65)


def supplement(output, model):
    fig, (box,) = s.canvas(165, 86, [161], [])
    s.panel_letter(fig, box, "a")
    s.panel_text(fig, box, 28, 2.3, "Nominal reference-rank decisions", va="center")
    for x, color, text in ((112, s.CROSS, "Crosslabel"), (136, s.PAIR, "Marked pair")):
        s.panel_text(fig, box, x, 2.3, text, color=color, va="center")
    for i, (x, label, limit) in enumerate(((31, "Trimer detection", 1.02),
                                           (77, "False alerts at references", .22),
                                           (123, "False alerts between references", .22))):
        ax = rate_axes(fig, box, x, 13, 34, 25, 4, limit,
                       [0, .25, .5, .75, 1] if i == 0 else [0, .05, .1, .2])
        s.panel_text(fig, box, x + 17, 10.4, label, ha="center", va="bottom")
        if i:
            ax.axvline(model["nominal"], color=s.REF, lw=.6, linestyle=":")
        for yy, case in enumerate(model["rates"]):
            if i == 0:
                ax.text(-.115, yy, case["label"], transform=ax.get_yaxis_transform(),
                        ha="right", va="center")
            if case["groups"][i] is not None:
                for row, color, mark, off in zip(case["groups"][i], (s.CROSS, s.PAIR),
                                                 ("o", "s"), (-.19, .19)):
                    rate_dot(ax, row, yy + off, color, mark)
    bottom = s.axis(fig, box, 31, 57, 66, 19)
    bottom.set(xlim=(0, 1.025), ylim=(2.55, -.55),
               xticks=[0, .25, .5, .75, 1.], yticks=[])
    bottom.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    bottom.spines["left"].set_visible(False)
    upper = box
    # The second row is its own GridSpec slot for the shared panel-letter offset.
    grid = fig.add_gridspec(1, 1, left=box.x0, right=box.x1, bottom=box.y0,
                           top=s.mm(fig, box, 0, 47)[1])
    lower = grid[0].get_position(fig)
    s.panel_letter(fig, lower, "b")
    s.panel_text(fig, upper, 28, 49.4, "Trimer detection at five false alerts per 100", va="center")
    from matplotlib.lines import Line2D
    for method, title, x, y, marker, fill in (
        ("crosslabel_lower", "Crosslabel lower endpoint", 109, 59, "o", s.CROSS),
        ("crosslabel_point", "Crosslabel point estimate", 109, 66, "s", "white"),
        ("dtc_connected", "Connected triple correlation", 109, 73, "^", s.TRIPLE),
    ):
        color, _ = METHODS[method]
        fig.add_artist(Line2D(*([v] for v in s.mm(fig, box, x, y)), transform=fig.transFigure,
                              marker=marker, color=color, mfc=fill, mew=.75,
                              markersize=3.1, linestyle="none"))
        s.panel_text(fig, box, x + 3, y, title, va="center")
    for yy, case in enumerate(model["matched"]):
        bottom.text(-.06, yy, case["label"], transform=bottom.get_yaxis_transform(),
                    ha="right", va="center")
        for method, marker, fill, dy in (
            ("crosslabel_lower", "o", s.CROSS, -.25),
            ("crosslabel_point", "s", "white", 0),
            ("dtc_connected", "^", s.TRIPLE, .25),
        ):
            bottom.plot(case["counts"][method] / model["matched_alternative_n"],
                        yy + dy, marker=marker, color=METHODS[method][0],
                        markersize=3.1, mfc=fill, mew=.75, linestyle="none")
    s.save(fig, output, "crosslabel_suppfig1")


def power(output, model):
    fig, boxes = s.canvas(165, 70, [76, 76], [9])
    colors = (s.PAIR, s.CROSS, s.TRIPLE)
    markers = ("o", "s", "^")
    handles = []
    for index, (box, alternative, label) in enumerate(zip(
        boxes, (True, False), ("Alternative detection (%)", "False alerts (%)")
    )):
        s.panel_letter(fig, box, "ab"[index])
        ax = s.axis(fig, box, 11, 6, 62, 45)
        rows = [r for r in model["cells"] if r["alternative"] == alternative]
        upper = 1.02 if alternative else max(.25, max(r["ci95"][1] for r in rows) * 1.05)
        ax.set(xlim=(.10, 1.04), ylim=(-.015, upper),
               xticks=model["retention_probabilities"], ylabel=label)
        ax.set_xlabel("Retention probability", labelpad=3)
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        if not alternative:
            ax.axhline(model["false_alert_budget"], color=s.REF, lw=.6, linestyle=":", zorder=0)
        for count, color, marker in zip(model["expected_record_counts"], colors, markers):
            group = sorted((r for r in rows if r["expected_records"] == count),
                           key=lambda r: r["retention"])
            if any(r["supported"] != r["planned"] for r in group):
                raise ValueError("Unsupported power cell needs a separately labelled display")
            x = np.array([r["retention"] for r in group])
            y = np.array([r["alerts"] / r["planned"] for r in group])
            bounds = np.array([r["ci95"] for r in group])
            ax.errorbar(x, y, yerr=np.array([y-bounds[:,0], bounds[:,1]-y]),
                        color=color, marker=marker, ms=3.1, lw=.8, elinewidth=.65,
                        capsize=1.3, capthick=.65, mfc="white")
            if index == 0:
                handles.append(Line2D([], [], color=color, marker=marker, mfc="white",
                                      label=f"{count:,} expected anchors"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .02),
               ncols=3, frameon=False, handlelength=2, columnspacing=1.5, borderaxespad=0)
    s.save(fig, output, "crosslabel_power")


def area(output, model):
    fig, boxes = s.canvas(165, 70, [76, 76], [9])
    handles = []
    for index, (box, alternative, label) in enumerate(zip(
        boxes, (True, False), ("Alternative detection (%)", "False alerts (%)")
    )):
        s.panel_letter(fig, box, "ab"[index])
        ax = s.axis(fig, box, 11, 6, 62, 45)
        rows = [r for r in model["groups"] if r["alternative"] == alternative]
        upper = 1.02 if alternative else max(.25, max(r["ci95"][1] for r in rows)*1.05)
        ax.set(xlim=(80, 420), ylim=(-.015, upper), xticks=model["areas_um2"], ylabel=label)
        ax.set_xlabel("Parent-window area (µm²)", labelpad=3)
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        if not alternative:
            ax.axhline(model["false_alert_budget"], color=s.REF, lw=.6, linestyle=":", zorder=0)
        for p, color, marker in zip(model["retention"], (s.PAIR, s.CROSS, s.TRIPLE), ("o", "s", "^")):
            group = sorted((r for r in rows if r["retention"] == p), key=lambda r: r["area"])
            if any(r["supported"] != r["planned"] for r in group):
                raise ValueError("Unsupported area-study cells need an explicit display")
            x = np.array([r["area"] for r in group])
            y = np.array([r["alerts"]/r["planned"] for r in group])
            ci = np.array([r["ci95"] for r in group])
            ax.errorbar(x, y, yerr=np.array([y-ci[:,0], ci[:,1]-y]),
                        color=color, marker=marker, ms=3.1, lw=.8, elinewidth=.65,
                        capsize=1.3, capthick=.65, mfc="white")
            if index == 0:
                handles.append(Line2D([], [], color=color, marker=marker, mfc="white",
                                      label=f"Retention {p:.2f}"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .02),
               ncols=3, frameon=False, handlelength=2, columnspacing=2, borderaxespad=0)
    s.save(fig, output, "crosslabel_area")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--compact-workflow", action="store_true")
    parser.add_argument("--coverage-in-supplement", action="store_true")
    parser.add_argument("--from-values", type=Path, help="Render an audited figure_values.json")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.from_values is None:
        import figure_data as data

        receipt = {
            "compact_workflow": args.compact_workflow,
            "coverage_in_supplement": args.coverage_in_supplement,
            "figure_1": data.workflow(),
            "figure_2": data.matched_scores(),
            "figure_3": data.roc_data(),
            "coverage": data.coverage_data(),
            "figure_4": data.nup96_data(),
            "supplement": data.supplement_rates(),
        }
        if (data.DATA / "retention-power-nov18/SUMMARY.json").exists():
            receipt["retention_power"] = data.retention_power()
        if (data.DATA / "area-power-nov18/SUMMARY.json").exists():
            receipt["area_power"] = data.area_power()
    else:
        if args.compact_workflow or args.coverage_in_supplement:
            parser.error("--from-values uses the layout recorded in the input")
        receipt = json.loads(args.from_values.read_text(encoding="utf-8"))
    workflow(args.output, receipt["figure_1"], receipt["compact_workflow"])
    matched(args.output, receipt["figure_2"])
    roc(args.output, receipt["figure_3"], receipt["coverage"], receipt["coverage_in_supplement"])
    nup96(args.output, receipt["figure_4"])
    supplement(args.output, receipt["supplement"])
    if "retention_power" in receipt:
        power(args.output, receipt["retention_power"])
    if "area_power" in receipt:
        area(args.output, receipt["area_power"])
    (args.output / "figure_values.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
