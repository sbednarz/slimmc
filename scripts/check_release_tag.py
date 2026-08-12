#!/usr/bin/env python3
"""Reject a component release tag that differs from its checked-in version."""
from __future__ import annotations

import argparse

from check_versions import check


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component", choices=("slimmc", "pyslimmc", "pyslimmc-opt"), required=True)
    parser.add_argument("--tag", required=True)
    args = parser.parse_args()

    current = check()
    expected = f"{args.component}-v{current[args.component]}"
    if args.tag != expected:
        raise SystemExit(
            f"Release tag mismatch: received {args.tag!r}, expected {expected!r}. "
            "Update the component version files or create the exact matching tag."
        )
    print(f"release tag: PASS ({expected})")


if __name__ == "__main__":
    main()
