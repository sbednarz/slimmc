#!/usr/bin/env python3
"""Update the build-time version copies for one independently released component."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace(relative: str, pattern: str, replacement: str) -> None:
    path = ROOT / relative
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"Expected exactly one version field in {relative}, found {count}")
    path.write_text(updated, encoding="utf-8")


def set_slimmc(version: str) -> None:
    for relative in ("VERSION", "homo/VERSION", "copo/VERSION"):
        (ROOT / relative).write_text(version + "\n", encoding="utf-8")
    for relative, pattern, prefix in (
        ("homo/src/slimmc_types.nim", r'(SlimmcVersion\*\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', r'\1'),
        ("copo/src/copo_types.nim", r'(SlimmcVersion\*\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', r'\1'),
        ("cli/slimmc_cli.nim", r'(CliVersion\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', r'\1'),
        ("cli/slimmc_summary.nim", r'(AppVersion\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', r'\1'),
        ("copo/slimmc-copo.nimble", r'^(version\s*=\s*)"[^"]+"', r'\1'),
        ("CITATION.cff", r'^(version:\s*)"[^"]+"', r'\1'),
    ):
        replace(relative, pattern, prefix + f'"{version}"')


def set_pyslimmc(version: str) -> None:
    replace("pyslimmc/_version.py", r'^(__version__\s*=\s*)"[^"]+"', rf'\1"{version}"')
    replace("pyproject.toml", r'^(version\s*=\s*)"[^"]+"', rf'\1"{version}"')
    replace("cli/slimmc_cli.nim", r'(PyslimmcVersion\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', rf'\1"{version}"')


def set_opt(version: str) -> None:
    replace("pyslimmc_opt/__init__.py", r'^(__version__\s*=\s*)"[^"]+"', rf'\1"{version}"')
    replace("pyslimmc_opt/pyproject.toml", r'^(version\s*=\s*)"[^"]+"', rf'\1"{version}"')
    replace("cli/slimmc_cli.nim", r'(PyslimmcOptVersion\s*\{\.strdefine\.\}\s*=\s*)"[^"]+"', rf'\1"{version}"')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("component", choices=("slimmc", "pyslimmc", "pyslimmc-opt"))
    parser.add_argument("version", help="Exact version stored in files and used after the tag's -v prefix")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9][0-9A-Za-z.-]*", args.version):
        raise SystemExit(f"Unsupported version syntax: {args.version!r}")
    {"slimmc": set_slimmc, "pyslimmc": set_pyslimmc, "pyslimmc-opt": set_opt}[args.component](args.version)
    print(f"updated {args.component} to {args.version}")
    print(f"expected tag: {args.component}-v{args.version}")


if __name__ == "__main__":
    main()
