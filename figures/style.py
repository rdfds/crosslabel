"""Print dimensions, typography and final-artwork checks."""
from pathlib import Path
import re
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
INK = "#222222"
REF = "#555555"
CROSS = "#0072B2"
TRIPLE = "#D55E00"
PAIR = "#7A7A7A"
ADEQUATE = "#B5B5B5"
ALTERNATIVE = "#222222"
ROUNDS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]


def configure():
    plt.style.use(HERE / "paper.mplstyle")
    mac_fonts = Path("/System/Library/Fonts/Supplemental")
    for name in ("Arial.ttf", "Arial Bold.ttf"):
        if (mac_fonts / name).exists():
            font_manager.fontManager.addfont(mac_fonts / name)
    for family in ("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans"):
        try:
            regular = font_manager.findfont(family, fallback_to_default=False)
            bold = font_manager.findfont(
                font_manager.FontProperties(family=family, weight="bold"),
                fallback_to_default=False,
            )
        except ValueError:
            continue
        mpl.rcParams.update({"font.family": family, "text.color": INK,
                             "axes.labelcolor": INK})
        return {font_manager.get_font(path).postscript_name for path in (regular, bold)}
    raise RuntimeError("Install Arial, Helvetica, Liberation Sans or Nimbus Sans.")


FONTS = configure()


def canvas(width, height, widths, gaps, *, bottom=3, top=3, left=2, right=2):
    """Construct panel slots in millimetres, including the explicit gaps."""
    if abs(sum(widths) + sum(gaps) - (width - left - right)) > 1e-8:
        raise ValueError("Panel sizes and margins must sum to the page width.")
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    columns = []
    for index, value in enumerate(widths):
        columns.append(value)
        if index < len(gaps):
            columns.append(gaps[index])
    grid = fig.add_gridspec(
        1, len(columns), left=left / width, right=1 - right / width,
        bottom=bottom / height, top=1 - top / height,
        width_ratios=columns, wspace=0,
    )
    return fig, [grid[0, index].get_position(fig) for index in range(0, len(columns), 2)]


def mm(fig, box, x, y):
    width, height = fig.get_size_inches() * 25.4
    return box.x0 + x / width, box.y1 - y / height


def panel_text(fig, box, x, y, label, **kwargs):
    return fig.text(*mm(fig, box, x, y), label, **kwargs)


def panel_letter(fig, box, label):
    panel_text(fig, box, 0.6, 0.4, label, fontsize=9, weight="bold", va="top")


def axis(fig, box, x, y, width, height, **kwargs):
    left, bottom = mm(fig, box, x, y + height)
    fw, fh = fig.get_size_inches() * 25.4
    return fig.add_axes([left, bottom, width / fw, height / fh], **kwargs)


def bar(ax, x, y, length, label):
    ax.plot([x, x + length], [y, y], color=INK, lw=1)
    ax.annotate(label, (x + length / 2, y), xytext=(0, 3),
                textcoords="offset points", ha="center", va="bottom")


def audit_text(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    page = fig.bbox
    texts = []
    for item in fig.findobj(Text):
        if not item.get_visible() or not item.get_text().strip():
            continue
        if item.get_fontsize() < 6.5:
            raise ValueError(f"Text below 6.5 pt: {item.get_text()!r}")
        bounds = item.get_window_extent(renderer)
        if min(bounds.x0 - page.x0, bounds.y0 - page.y0,
               page.x1 - bounds.x1, page.y1 - bounds.y1) < -0.3:
            raise ValueError(f"Text outside figure: {item.get_text()!r}")
        for other, extent in texts:
            # A subpixel touch is not an overlap in print.
            overlap_x = min(bounds.x1, extent.x1) - max(bounds.x0, extent.x0)
            overlap_y = min(bounds.y1, extent.y1) - max(bounds.y0, extent.y0)
            if overlap_x > 0.6 and overlap_y > 0.6:
                raise ValueError(f"Text overlap: {other!r} and {item.get_text()!r}")
        texts.append((item.get_text(), bounds))


def save(fig, output, name):
    audit_text(fig)
    pdf = output / f"{name}.pdf"
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(output / f"{name}.png", dpi=300)
    result = subprocess.check_output(["pdffonts", str(pdf)], text=True)
    found = {re.sub(r"^[A-Z]{6}\+", "", row.split()[0])
             for row in result.splitlines()[2:] if row.strip()}
    if not found or not found <= FONTS:
        raise ValueError(f"{name}: expected {FONTS}, found {found}")
    plt.close(fig)
