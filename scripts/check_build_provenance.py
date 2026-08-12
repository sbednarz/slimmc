#!/usr/bin/env python3
"""Verify that a built Slimmc executable exposes compile-time provenance."""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("executable", type=Path)
    ap.add_argument("--require-git", action="store_true")
    args = ap.parse_args()

    proc = subprocess.run([str(args.executable), "--version"], capture_output=True, text=True, check=True)
    text = proc.stdout + proc.stderr
    required = ("build timestamp:",)
    missing = [label for label in required if label not in text.lower()]
    if missing:
        raise SystemExit(f"{args.executable}: missing build provenance fields: {', '.join(missing)}")
    if args.require_git:
        m = re.search(r"git commit:\s*(\S+)", text, flags=re.I)
        if not m or not re.fullmatch(r"[0-9a-fA-F]{7,64}", m.group(1)):
            raise SystemExit(f"{args.executable}: release build has no valid compiled-in Git commit SHA")
        if not re.search(r"git dirty:\s*false", text, flags=re.I):
            raise SystemExit(f"{args.executable}: release build is not from a clean Git tree")
    print("build provenance: PASS")


if __name__ == "__main__":
    main()
