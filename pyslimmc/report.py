"""Small Matplotlib-only multipage PDF report builder.

The module intentionally implements a flow of only a few content blocks.  It
is not a document-layout engine.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
from pathlib import Path
import textwrap
from typing import Any, Callable


_ALIGNMENTS = {"left", "center", "right"}
_MATH_FONTS = {"dejavusans", "dejavuserif", "cm", "stix", "stixsans"}
_A4 = (8.27, 11.69)
_PLOT_TOP = 0.28
_PLOT_BOTTOM = 0.34


@dataclass(frozen=True)
class _Block:
    kind: str
    value: Any
    options: dict[str, Any]


class Report:
    """A simple flowing A4 report rendered entirely by Matplotlib.

    Methods append blocks and return ``self``, so calls may be chained.  Plot
    blocks accept a ``callable(ax)``, a pyslimmc object exposing
    ``plot(ax=...)``, or a Matplotlib Figure.  Figure-only plot APIs are
    rasterized into their report block; ``plot(ax=...)`` remains vector.
    """

    def __init__(self, title: str | None = None, *, orientation: str = "portrait",
                 font: str = "DejaVu Sans", font_size: float = 10.0,
                 math_font: str = "dejavusans"):
        if orientation not in {"portrait", "landscape"}:
            raise ValueError("orientation must be 'portrait' or 'landscape'")
        if not isinstance(font, str) or not font.strip():
            raise ValueError("font must be a non-empty family name")
        if isinstance(font_size, bool) or float(font_size) <= 0:
            raise ValueError("font_size must be positive")
        self.title = None if title is None else str(title)
        self.orientation = orientation
        self.font = font
        self.font_size = float(font_size)
        self.math_font = _math_font(math_font)
        self._blocks: list[_Block] = []

    def text(self, value: Any, *, size: float | None = None, font: str | None = None,
             align: str = "left", weight: str = "normal") -> "Report":
        """Append wrapping plain text."""
        size = self.font_size if size is None else _positive(size, "size")
        self._blocks.append(_Block("text", str(value), {
            "size": size, "font": font or self.font,
            "align": _alignment(align), "weight": str(weight),
        }))
        return self

    def text_raw(self, value: Any, *, size: float | None = None,
                 font: str = "DejaVu Sans Mono", align: str = "left",
                 weight: str = "normal") -> "Report":
        """Append preformatted text without wrapping or whitespace changes.

        Newlines and leading spaces are preserved. Tabs are expanded to four
        spaces for deterministic rendering. Long lines may extend past the
        right margin; that is the intentional meaning of ``raw``.
        """
        size = self.font_size if size is None else _positive(size, "size")
        if not isinstance(font, str) or not font.strip():
            raise ValueError("font must be a non-empty family name")
        self._blocks.append(_Block("text_raw", str(value).expandtabs(4), {
            "size": size, "font": font, "align": _alignment(align),
            "weight": str(weight),
        }))
        return self

    def math(self, expression: str, *, size: float | None = None,
             align: str = "center", font: str | None = None) -> "Report":
        """Append one Matplotlib MathText expression (no external LaTeX)."""
        if not isinstance(expression, str) or not expression.strip():
            raise ValueError("expression must be a non-empty string")
        value = expression.strip()
        if not (value.startswith("$") and value.endswith("$")):
            value = f"${value}$"
        size = self.font_size * 1.25 if size is None else _positive(size, "size")
        self._blocks.append(_Block("math", value, {
            "size": size, "align": _alignment(align),
            "font": self.math_font if font is None else _math_font(font),
        }))
        return self

    def vspace(self, lines: float = 1.0) -> "Report":
        """Append vertical whitespace measured in default text lines."""
        lines = _positive(lines, "lines")
        self._blocks.append(_Block("vspace", None, {"lines": lines}))
        return self

    def plot(self, value: Any, *, height: float | None = None,
             span: str | None = None, align: str = "center", **kwargs: Any) -> "Report":
        """Append a plot block.

        ``value`` may be ``callable(ax)``, an object with ``plot()``, or an
        existing Matplotlib Figure.  ``kwargs`` are forwarded to the callable
        or object's ``plot`` method.
        """
        if not callable(value) and not hasattr(value, "plot") and not _is_figure(value):
            raise TypeError("plot value must be callable, expose plot(), or be a Matplotlib Figure")
        align = _alignment(align)
        if span is None:
            height = 3.2 if height is None else _positive(height, "height")
            options = {"height": height, "width": None, "span": None,
                       "align": align, "kwargs": dict(kwargs)}
        else:
            if height is not None:
                raise ValueError("height cannot be combined with span; the publication preset fixes both dimensions")
            from .plotting import figure_size
            style = kwargs.get("style", "screen")
            width, fixed_height = figure_size(style, span=span)
            options = {"height": fixed_height, "width": width, "span": span,
                       "align": align, "kwargs": dict(kwargs)}
        self._blocks.append(_Block("plot", value, options))
        return self

    def page_break(self) -> "Report":
        """Finish the current page.  A trailing break creates no blank page."""
        self._blocks.append(_Block("page_break", None, {}))
        return self

    def save(self, path: str | Path) -> Path:
        """Render the report to an atomically replaced multipage PDF."""
        try:
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
        except ImportError as exc:
            raise ImportError("Report.save() requires matplotlib") from exc

        target = Path(path)
        if target.suffix.lower() != ".pdf":
            raise ValueError("report path must end with .pdf")
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(target.name + ".tmp")
        page_size = _A4 if self.orientation == "portrait" else _A4[::-1]
        state: dict[str, Any] = {}
        page_number = 0

        def new_page() -> None:
            nonlocal page_number
            page_number += 1
            figure = plt.figure(figsize=page_size)
            state.update(figure=figure, cursor=page_size[1] - 0.72, used=False)
            if self.title:
                figure.text(0.08, 0.965, self.title, ha="left", va="top",
                            family=self.font, fontsize=self.font_size * 1.4,
                            fontweight="bold")
                state["cursor"] -= 0.42

        def finish_page(pdf: Any) -> None:
            figure = state["figure"]
            figure.text(0.92, 0.025, str(page_number), ha="right", va="bottom",
                        family=self.font, fontsize=max(7.0, self.font_size * 0.8), color="#555555")
            pdf.savefig(figure)
            plt.close(figure)

        try:
            with PdfPages(temporary) as pdf:
                new_page()
                for block in self._blocks:
                    if block.kind == "page_break":
                        if state["used"]:
                            finish_page(pdf)
                            new_page()
                        continue
                    prepared, height = self._prepare(block, page_size)
                    bottom = 0.62
                    if state["used"] and state["cursor"] - height < bottom:
                        finish_page(pdf)
                        new_page()
                    if height > state["cursor"] - bottom:
                        raise ValueError(f"{block.kind} block is taller than one report page")
                    self._render(state["figure"], block, prepared, state["cursor"], height, page_size)
                    state["cursor"] -= height
                    state["used"] = True
                # A new page after a trailing page_break is intentionally not saved.
                if state["used"] or page_number == 1:
                    finish_page(pdf)
                else:
                    plt.close(state["figure"])
            temporary.replace(target)
        except Exception:
            if "figure" in state:
                plt.close(state["figure"])
            temporary.unlink(missing_ok=True)
            raise
        return target

    def _prepare(self, block: _Block, page_size: tuple[float, float]) -> tuple[Any, float]:
        if block.kind == "plot":
            if block.options["span"] is not None:
                return None, block.options["height"] + 0.12
            return None, block.options["height"] + _PLOT_TOP + _PLOT_BOTTOM
        if block.kind == "vspace":
            return None, block.options["lines"] * self.font_size / 72.0 * 1.35
        size = block.options["size"]
        if block.kind == "math":
            return block.value, size / 72.0 * 1.9 + 0.10
        if block.kind == "text_raw":
            lines = block.value.splitlines() or [""]
            return block.value, len(lines) * size / 72.0 * 1.35 + 0.10
        usable_width = page_size[0] * 0.84
        chars = max(12, int(usable_width * 72.0 / (size * 0.53)))
        lines: list[str] = []
        for paragraph in block.value.splitlines() or [""]:
            lines.extend(textwrap.wrap(paragraph, width=chars, replace_whitespace=False,
                                       drop_whitespace=True) or [""])
        height = len(lines) * size / 72.0 * 1.35 + 0.10
        return "\n".join(lines), height

    def _render(self, figure: Any, block: _Block, prepared: Any, cursor: float,
                height: float, page_size: tuple[float, float]) -> None:
        y_top = cursor / page_size[1]
        align = block.options.get("align", "left")
        x = {"left": 0.08, "center": 0.5, "right": 0.92}[align]
        if block.kind in {"text", "text_raw"}:
            figure.text(x, y_top, prepared, ha=align, va="top",
                        family=block.options["font"], fontsize=block.options["size"],
                        fontweight=block.options["weight"], linespacing=1.35)
        elif block.kind == "math":
            figure.text(x, y_top - height / page_size[1] * 0.18, prepared,
                        ha=align, va="top", fontsize=block.options["size"],
                        math_fontfamily=block.options["font"])
        elif block.kind == "plot":
            plot_height = block.options["height"]
            plot_width = block.options["width"]
            if plot_width is not None:
                from .plotting import figure_margins
                left_margin, right_margin, bottom_margin, top_margin = figure_margins()
                usable_left = page_size[0] * 0.08
                usable_right = page_size[0] * 0.92
                align = block.options["align"]
                canvas_left = {
                    "left": usable_left,
                    "center": (page_size[0] - plot_width) / 2.0,
                    "right": usable_right - plot_width,
                }[align]
                canvas_bottom = cursor - plot_height
                ax = figure.add_axes((
                    (canvas_left + left_margin) / page_size[0],
                    (canvas_bottom + bottom_margin) / page_size[1],
                    (plot_width - left_margin - right_margin) / page_size[0],
                    (plot_height - bottom_margin - top_margin) / page_size[1],
                ))
                self._draw_plot(ax, block.value, block.options["kwargs"])
                return
            # Reserve a small cap above the Axes: Matplotlib draws an Axes title
            # outside its rectangle, so placing the rectangle directly at the
            # flow cursor would let the title touch the preceding block.
            bottom = (cursor - _PLOT_TOP - plot_height) / page_size[1]
            ax = figure.add_axes((0.10, bottom, 0.82, plot_height / page_size[1]))
            self._draw_plot(ax, block.value, block.options["kwargs"])

    @staticmethod
    def _draw_plot(ax: Any, value: Any, kwargs: dict[str, Any]) -> None:
        import matplotlib.pyplot as plt
        if _is_figure(value):
            _figure_image(ax, value)
            return
        if callable(value) and not hasattr(value, "plot"):
            result = value(ax, **kwargs)
        else:
            method: Callable[..., Any] = value.plot
            try:
                accepts_ax = "ax" in inspect.signature(method).parameters
            except (TypeError, ValueError):
                accepts_ax = True
            result = method(ax=ax, **kwargs) if accepts_ax else method(**kwargs)
        if _is_figure(result) and result is not ax.figure:
            _figure_image(ax, result)
            plt.close(result)


def _figure_image(ax: Any, figure: Any) -> None:
    """Place a Figure-only plot in a report block using Matplotlib's RGBA buffer."""
    figure.canvas.draw()
    image = figure.canvas.buffer_rgba()
    ax.imshow(image)
    ax.set_axis_off()


def _is_figure(value: Any) -> bool:
    try:
        from matplotlib.figure import Figure
    except ImportError:
        return False
    return isinstance(value, Figure)


def _positive(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a positive number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a positive number") from exc
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _alignment(value: str) -> str:
    if value not in _ALIGNMENTS:
        raise ValueError("align must be 'left', 'center', or 'right'")
    return value


def _math_font(value: str) -> str:
    if value not in _MATH_FONTS:
        raise ValueError("math font must be one of: " + ", ".join(sorted(_MATH_FONTS)))
    return value


def report(*args: Any, **kwargs: Any) -> Report:
    """Convenience constructor equivalent to ``Report(...)``."""
    return Report(*args, **kwargs)
