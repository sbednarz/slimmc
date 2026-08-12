from __future__ import annotations

import argparse
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--platform", required=True)
    parser.add_argument("--exe", default="")
    args = parser.parse_args()

    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    name = f"slimmc-{version}-{args.platform}"
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    target = dist / f"{name}.zip"
    files = [
        (ROOT / f"bin/slimmc{args.exe}", f"{name}/bin/slimmc{args.exe}"),
        (ROOT / f"bin/slimmc-summary{args.exe}", f"{name}/bin/slimmc-summary{args.exe}"),
        (ROOT / "README.md", f"{name}/README.md"),
        (ROOT / "LICENSE", f"{name}/LICENSE"),
        (ROOT / "VERSION", f"{name}/VERSION"),
    ]
    missing = [str(path) for path, _ in files if not path.is_file()]
    if missing:
        raise SystemExit("Missing release inputs: " + ", ".join(missing))
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for source, arcname in files:
            archive.write(source, arcname)
    print(target)


if __name__ == "__main__":
    main()
