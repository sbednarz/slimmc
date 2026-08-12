#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "docs" / "reference"
TARGETS = {
    "pyslimmc": REF / "PYSLIMMC_SIGNATURES.md",
    "pyslimmc_opt": REF / "PYSLIMMC_OPT_SIGNATURES.md",
}
sys.path.insert(0, str(ROOT))

PYSLIMMC_MODULES = (
    "pyslimmc.run",
    "pyslimmc.runs",
    "pyslimmc._storage",
    "pyslimmc.chains",
    "pyslimmc.chain_counts",
    "pyslimmc.spectra",
    "pyslimmc.composition_analysis",
    "pyslimmc.copolymerization",
    "pyslimmc.storage_analysis",
    "pyslimmc.report",
    "pyslimmc.summary",
    "pyslimmc.table",
    "pyslimmc.plotting",
)
PYSLIMMC_OPT_MODULES = ("pyslimmc_opt.study",)

PYSLIMMC_PREAMBLE = """# pyslimmc callable signatures

This is the exhaustive callable inventory for objects intentionally returned
to users by `pyslimmc`. It is generated from the installed source by
`scripts/update_api_signatures.py`; CI rejects a stale inventory. The tutorial
and semantic reference remain in [`../PYSLIMMC.md`](../PYSLIMMC.md) and
[`PYSLIMMC_API.md`](PYSLIMMC_API.md).

Classes exported at package root are also typing/inspection names. Construct
runs and collections with `open()` and `scan()`; do not instantiate their
Storage implementation classes directly.

## Shared parameter contract

| Parameter | Accepted values and meaning |
|---|---|
| `snapshot` | `"final"` (default), `"last"`, an integer snapshot ID, or a `StorageSnapshot`. |
| `pool` | `"all"`, `"live"`, `"dead"`, or a kinetic pool name; a sequence requests a grouped result where supported. |
| `series` | Mapping/name-to-population selectors for several distributions; mutually exclusive with a non-`all` `pool`. |
| `mass_model` | `"repeat_units"` or `"with_end_groups"`; `None` uses stored/default mass semantics. |
| `progress` | `None` uses `pyslimmc.options.progress`; `True` forces and `False` suppresses progress. |
| `method` | Distribution representation: `"sticks"`, `"hist"`, `"gaussian"`, or `"kde"`. |
| `basis` | `"number"` or `"mass"`. |
| `coordinate` | `"linear"` or `"log10"`. |
| `output` | `"amount"`, `"fraction"`, or `"density"`. |
| `normalization` | `"absolute"`, `"per_series"`, `"combined"`, or `"reference"`. |
| `bins`, `bin_width` | Alternative grid controls; they are mutually exclusive and must be positive. |
| `sigma` | Positive smoothing width for Gaussian/KDE methods; units follow `coordinate`. |
| `grid_step` | Positive output-grid step; units follow `coordinate`. |
| `reference` | Reference series name when `normalization="reference"`. |
| `style` | Plot style from `available_styles()`; default `"screen"`. |
| `span` | `None`, `"column"`, or `"double"`; controls owned figure geometry. |
| `ax` | Existing Matplotlib axes. If supplied, `span` must be omitted. |
| `path` | Optional output path. Methods return the written `Path` where applicable. |
| `dpi` | Positive raster resolution, normally `300`. |
| `metadata` | Export metadata mode; `"comments"` writes a commented header. |
| `layout` | Export layout `"wide"` or `"long"` where supported. |

For MWD the defaults are `basis="mass"`, `method="gaussian"`,
`coordinate="log10"`, `output="density"`, and
`normalization="per_series"`. For CLD they are `basis="number"`,
`method="sticks"`, `coordinate="linear"`, `output="fraction"`, and
`normalization="per_series"`. A neutral chain-mass spectrum uses exact sticks;
`normalize` is `"count"`, `"fraction"`, or `"base_peak"`.

## Top-level callables
"""

PYSLIMMC_OPT_PREAMBLE = """# pyslimmc-opt callable signatures

This is the exhaustive callable inventory for the public `pyslimmc_opt` API.
It is generated from the installed source by `scripts/update_api_signatures.py`;
CI rejects a stale inventory. The human-readable reference remains in
[`PYSLIMMC_OPT.md`](PYSLIMMC_OPT.md).

## Top-level callables
"""


def public_classes(module_name: str) -> list[type]:
    module = importlib.import_module(module_name)
    classes = []
    for name, value in vars(module).items():
        if name.startswith("_") or not inspect.isclass(value):
            continue
        if value.__module__ == module_name:
            classes.append(value)
    return sorted(classes, key=lambda item: item.__name__)


def callable_signature(value, *, drop_bound: bool = False) -> str:
    try:
        result = inspect.signature(value)
        parameters = tuple(result.parameters.values())
        if drop_bound and parameters and parameters[0].name in {"self", "cls"}:
            result = result.replace(parameters=parameters[1:])
        return str(result)
    except (TypeError, ValueError):
        return "(signature unavailable)"


def build_package(package_name: str, modules: tuple[str, ...], preamble: str) -> str:
    package = importlib.import_module(package_name)
    lines = [preamble.rstrip()]
    for name in package.__all__:
        value = getattr(package, name)
        if inspect.isfunction(value) or inspect.isclass(value):
            lines.append(f"- `{package.__name__}.{name}{callable_signature(value)}`")

    lines.extend(["", "## Object callables", ""])
    for module_name in modules:
        for cls in public_classes(module_name):
            methods = []
            properties = []
            for name, raw in vars(cls).items():
                if name.startswith("_"):
                    continue
                if isinstance(raw, property):
                    properties.append(name)
                elif callable(raw):
                    methods.append((name, callable_signature(raw, drop_bound=True)))
            if not methods and not properties:
                continue
            lines.append(f"### `{cls.__name__}`")
            lines.append("")
            if properties:
                lines.append("Properties: " + ", ".join(f"`{name}`" for name in sorted(properties)) + ".")
                lines.append("")
            for name, sig in sorted(methods):
                lines.append(f"- `{cls.__name__}.{name}{sig}`")
            lines.append("")
    if package_name == "pyslimmc":
        lines.extend([
            "",
            "## See also",
            "",
            "- [`../PYSLIMMC.md`](../PYSLIMMC.md) — task-oriented analysis guide",
            "- [`PYSLIMMC_API.md`](PYSLIMMC_API.md) — semantic API reference",
        ])
    else:
        lines.extend([
            "",
            "## See also",
            "",
            "- [`PYSLIMMC_OPT.md`](PYSLIMMC_OPT.md) — pyslimmc-opt guide and API semantics",
            "- [`../PYSLIMMC.md`](../PYSLIMMC.md) — analysis of Slimmc run results",
        ])
    return "\n".join(lines).rstrip() + "\n"


def builds() -> dict[str, str]:
    return {
        "pyslimmc": build_package("pyslimmc", PYSLIMMC_MODULES, PYSLIMMC_PREAMBLE),
        "pyslimmc_opt": build_package("pyslimmc_opt", PYSLIMMC_OPT_MODULES, PYSLIMMC_OPT_PREAMBLE),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = builds()
    REF.mkdir(parents=True, exist_ok=True)
    if args.check:
        stale = []
        for name, target in TARGETS.items():
            actual = target.read_text(encoding="utf-8") if target.is_file() else ""
            if actual != expected[name]:
                stale.append(str(target.relative_to(ROOT)))
        if stale:
            for path in stale:
                print(f"stale generated API signature reference: {path}")
            return 1
        print("API signature inventories: PASS")
        return 0
    for name, target in TARGETS.items():
        target.write_text(expected[name], encoding="utf-8")
        print(target.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
