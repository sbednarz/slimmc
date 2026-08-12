#!/usr/bin/env python3
"""Verify the canonical version and every build-time copy in the monorepo."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _text(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _line(relative: str) -> str:
    return _text(relative).strip()


def _match(relative: str, pattern: str) -> str:
    match = re.search(pattern, _text(relative), re.MULTILINE)
    if not match:
        raise SystemExit(f"Cannot read version from {relative}")
    return match.group(1)


def versions() -> dict[str, str]:
    return {
        "slimmc": _line("VERSION"),
        "pyslimmc": _match("pyslimmc/_version.py", r'^__version__\s*=\s*"([^"]+)"'),
        "pyslimmc-opt": _match("pyslimmc_opt/__init__.py", r'^__version__\s*=\s*"([^"]+)"'),
    }


def _require(label: str, actual: str, expected: str) -> None:
    if actual != expected:
        raise SystemExit(f"Version mismatch: {label} has {actual!r}, expected {expected!r}")


def check() -> dict[str, str]:
    current = versions()
    slimmc = current["slimmc"]
    pyslimmc = current["pyslimmc"]
    opt = current["pyslimmc-opt"]

    for path in ("homo/VERSION", "copo/VERSION"):
        _require(path, _line(path), slimmc)
    for path, pattern in (
        ("homo/src/slimmc_types.nim", r'SlimmcVersion\*\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        ("copo/src/copo_types.nim", r'SlimmcVersion\*\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        ("cli/slimmc_cli.nim", r'CliVersion\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        ("cli/slimmc_summary.nim", r'AppVersion\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        ("copo/slimmc-copo.nimble", r'^version\s*=\s*"([^"]+)"'),
        ("CITATION.cff", r'^version:\s*"([^"]+)"'),
    ):
        _require(path, _match(path, pattern), slimmc)

    _require("pyproject.toml", _match("pyproject.toml", r'^version\s*=\s*"([^"]+)"'), pyslimmc)
    _require(
        "cli/slimmc_cli.nim PyslimmcVersion",
        _match("cli/slimmc_cli.nim", r'PyslimmcVersion\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        pyslimmc,
    )
    _require(
        "pyslimmc_opt/pyproject.toml",
        _match("pyslimmc_opt/pyproject.toml", r'^version\s*=\s*"([^"]+)"'),
        opt,
    )
    _require(
        "cli/slimmc_cli.nim PyslimmcOptVersion",
        _match("cli/slimmc_cli.nim", r'PyslimmcOptVersion\s*\{\.strdefine\.\}\s*=\s*"([^"]+)"'),
        opt,
    )
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    current = check()
    if not args.quiet:
        for component, version in current.items():
            print(f"{component}={version}")
        print("version files: PASS")


if __name__ == "__main__":
    main()
