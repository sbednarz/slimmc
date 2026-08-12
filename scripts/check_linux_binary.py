"""Verify the portability contract of Linux release executables."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


GLIBC_RE = re.compile(r"GLIBC_(\d+)\.(\d+)")


def readelf(*args: str, binary: Path) -> str:
    completed = subprocess.run(
        ["readelf", "-W", *args, str(binary)],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout


def check_glibc(binary: Path, ceiling: tuple[int, int]) -> None:
    version_info = readelf("--version-info", binary=binary)
    versions = {tuple(map(int, match)) for match in GLIBC_RE.findall(version_info)}
    if not versions:
        raise SystemExit(f"{binary}: no GLIBC symbol versions found")
    required = max(versions)
    if required > ceiling:
        raise SystemExit(
            f"{binary}: requires GLIBC_{required[0]}.{required[1]}, "
            f"above GLIBC_{ceiling[0]}.{ceiling[1]}"
        )
    print(f"{binary}: maximum required ABI is GLIBC_{required[0]}.{required[1]}")


def check_static(binary: Path) -> None:
    program_headers = readelf("--program-headers", binary=binary)
    dynamic = readelf("--dynamic", binary=binary)
    if "INTERP" in program_headers:
        raise SystemExit(f"{binary}: contains a dynamic-program interpreter")
    if re.search(r"\(NEEDED\)", dynamic):
        raise SystemExit(f"{binary}: contains a dynamically needed library")
    print(f"{binary}: static ELF (no INTERP and no NEEDED entries)")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--glibc-max", metavar="MAJOR.MINOR")
    mode.add_argument("--static", action="store_true")
    parser.add_argument("binaries", nargs="+", type=Path)
    args = parser.parse_args()

    for binary in args.binaries:
        if not binary.is_file():
            raise SystemExit(f"missing binary: {binary}")
        if args.static:
            check_static(binary)
        else:
            ceiling = tuple(map(int, args.glibc_max.split(".")))
            if len(ceiling) != 2:
                raise SystemExit("--glibc-max must have the form MAJOR.MINOR")
            check_glibc(binary, ceiling)


if __name__ == "__main__":
    main()
