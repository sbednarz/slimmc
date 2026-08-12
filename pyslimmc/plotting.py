"""Local plotting styles used by distributions and plots.

The module deliberately never mutates :mod:`matplotlib` ``rcParams``.  A
style is applied only to the figure, axes and artists created by one call.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any


@dataclass(frozen=True)
class PlotStyle:
    name: str
    palette: tuple[str, ...]
    foreground: str
    background: str
    axes_background: str
    grid_color: str
    font_family: str
    font_size: float
    title_size: float
    line_width: float
    grid: bool
    column_size: tuple[float, float]
    double_size: tuple[float, float]
    default_span: str


_MM = 1.0 / 25.4
_COLUMN_SIZE = (82.0 * _MM, 62.0 * _MM)
_DOUBLE_SIZE = (170.0 * _MM, 90.0 * _MM)
_FIGURE_MARGINS = (0.62, 0.12, 0.50, 0.28)  # left, right, bottom, top / inch


_STYLES = MappingProxyType({
    "publication_bw": PlotStyle(
        "publication_bw", ("#000000", "#555555", "#888888", "#bbbbbb"),
        "#000000", "#ffffff", "#ffffff", "#c8c8c8", "DejaVu Sans",
        8.0, 9.0, 1.15, False, _COLUMN_SIZE, _DOUBLE_SIZE, "column",
    ),
    "publication_color": PlotStyle(
        "publication_color",
        ("#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9"),
        "#111111", "#ffffff", "#ffffff", "#d9d9d9", "DejaVu Sans",
        8.0, 9.0, 1.25, True, _COLUMN_SIZE, _DOUBLE_SIZE, "column",
    ),
    "screen": PlotStyle(
        "screen",
        ("#3366cc", "#dc3912", "#109618", "#990099", "#ff9900", "#0099c6"),
        "#202124", "#ffffff", "#ffffff", "#dfe3e8", "DejaVu Sans",
        10.0, 11.0, 1.5, True, _COLUMN_SIZE, _DOUBLE_SIZE, "double",
    ),
})


def available_styles() -> tuple[str, ...]:
    return tuple(_STYLES)


def get_style(name: str = "screen") -> PlotStyle:
    try:
        return _STYLES[name]
    except KeyError as exc:
        raise ValueError(
            f"unknown plot style {name!r}; available: {', '.join(available_styles())}"
        ) from exc


def style_kwargs(name: str = "screen", *, index: int = 0) -> dict[str, Any]:
    style = get_style(name)
    return {
        "color": style.palette[index % len(style.palette)],
        "linewidth": style.line_width,
    }


def figure_size(name: str = "screen", *, span: str | None = None) -> tuple[float, float]:
    """Return the exact physical figure size in inches.

    ``publication_*`` defaults to an A4 single-column 82 x 62 mm figure;
    ``screen`` defaults to the 170 x 90 mm double-column canvas.  Passing
    ``span=`` always wins over the style default.
    """
    style = get_style(name)
    selected = style.default_span if span is None else span
    if selected == "column":
        return style.column_size
    if selected == "double":
        return style.double_size
    raise ValueError("span must be 'column' or 'double'")


def create_axes(name: str = "screen", *, span: str | None = None):
    """Create a Figure/Axes pair with fixed publication geometry.

    Margins are physical rather than proportional, so labels occupy the same
    amount of paper in a one- and two-column figure.  No automatic layout
    engine or global rcParams mutation is involved.
    """
    import matplotlib.pyplot as plt
    width, height = figure_size(name, span=span)
    left, right, bottom, top = figure_margins()
    figure, ax = plt.subplots(figsize=(width, height))
    figure.subplots_adjust(
        left=left / width,
        right=1.0 - right / width,
        bottom=bottom / height,
        top=1.0 - top / height,
    )
    return figure, ax


def figure_margins() -> tuple[float, float, float, float]:
    """Physical left/right/bottom/top margins used by standalone figures."""
    return _FIGURE_MARGINS


def require_owned_geometry(ax, span: str | None) -> None:
    """Reject a misleading span when the caller owns the supplied Axes."""
    if ax is not None and span is not None:
        raise ValueError("span controls figures created by plot(); omit span when ax is supplied")


def apply_axes_style(ax, name: str = "screen"):
    """Apply one style locally and return ``ax``.

    Existing labels and artists are retained.  No global Matplotlib state is
    read or written beyond normal figure construction.
    """
    style = get_style(name)
    figure = ax.figure
    figure.set_facecolor(style.background)
    ax.set_facecolor(style.axes_background)
    ax.tick_params(colors=style.foreground, labelsize=style.font_size)
    for spine in ax.spines.values():
        spine.set_color(style.foreground)
        spine.set_linewidth(0.75)
    ax.xaxis.label.set_color(style.foreground)
    ax.yaxis.label.set_color(style.foreground)
    ax.title.set_color(style.foreground)
    for text in (ax.xaxis.label, ax.yaxis.label, ax.title, *ax.get_xticklabels(), *ax.get_yticklabels()):
        text.set_fontfamily(style.font_family)
    ax.xaxis.label.set_fontsize(style.font_size)
    ax.yaxis.label.set_fontsize(style.font_size)
    ax.title.set_fontsize(style.title_size)
    if style.grid:
        ax.grid(True, color=style.grid_color, linewidth=0.6, alpha=0.75)
    else:
        ax.grid(False)
    legend = ax.get_legend()
    if legend is not None:
        legend.get_frame().set_edgecolor(style.grid_color)
        legend.get_frame().set_facecolor(style.axes_background)
        for text in legend.get_texts():
            text.set_color(style.foreground)
            text.set_fontfamily(style.font_family)
            text.set_fontsize(style.font_size)
    return ax
