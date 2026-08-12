#!/usr/bin/env python
"""Small, dependency-free Python environment preflight for Make targets."""
from __future__ import annotations

import argparse
import importlib.util
import os
import platform
import sys
from pathlib import Path


def module_location(name: str) -> str:
    spec = importlib.util.find_spec(name)
    if spec is None:
        return "MISSING"
    if spec.origin:
        return spec.origin
    if spec.submodule_search_locations:
        return os.pathsep.join(str(Path(p)) for p in spec.submodule_search_locations)
    return "available"


def main() -> int:
    parser = argparse.ArgumentParser(description="Check the Python interpreter used by Slimmc")
    parser.add_argument("--require", nargs="*", default=[], metavar="MODULE")
    args = parser.parse_args()

    print(f"Python executable: {sys.executable}")
    print(f"Python version:    {platform.python_version()}")
    print(f"Virtual env:       {sys.prefix}")
    print(f"PYTHONPATH:        {os.environ.get('PYTHONPATH', '<unset>')}")

    missing: list[str] = []
    for name in args.require:
        location = module_location(name)
        print(f"{name:<18} {location}")
        if location == "MISSING":
            missing.append(name)

    if missing:
        joined = " ".join(missing)
        print("", file=sys.stderr)
        print("Missing Python modules: " + ", ".join(missing), file=sys.stderr)
        print(f"Install test dependencies with this exact interpreter:\n  {sys.executable} -m pip install -e '.[test]'", file=sys.stderr)
        print(f"Or install only the missing modules:\n  {sys.executable} -m pip install {joined}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
