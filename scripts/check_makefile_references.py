#!/usr/bin/env python3
"""Reject literal source/script paths in Makefiles that do not exist."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUFFIXES = (".py", ".nim", ".model")


def literal_references(path: Path) -> set[str]:
    references: set[str] = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0]
        for token in re.split(r"[\s:=]+", line):
            token = token.strip("'\"(),;\\")
            if not token.endswith(SUFFIXES):
                continue
            if "$" in token or "*" in token or token.startswith("-"):
                continue
            references.add(token)
    return references


def main() -> int:
    errors: list[str] = []
    makefiles = sorted(ROOT.rglob("Makefile"))
    for makefile in makefiles:
        for reference in sorted(literal_references(makefile)):
            resolved = (makefile.parent / reference).resolve()
            if not resolved.exists():
                errors.append(
                    f"{makefile.relative_to(ROOT)} references missing {reference}"
                )
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Makefile references: {len(makefiles)} file(s) PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
