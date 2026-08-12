from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FAMILY_ROOT = ROOT.parent


def remove_results(base: Path) -> None:
    if not base.exists():
        return
    for path in sorted(base.rglob("results*"), reverse=True):
        if path.is_dir():
            print("removing", path)
            shutil.rmtree(path)


def clear_out_dirs(base: Path) -> None:
    if not base.exists():
        return
    for path in base.rglob("out"):
        if not path.is_dir():
            continue
        for child in list(path.iterdir()):
            if child.name in {".keep", ".gitkeep"}:
                continue
            print("removing", child)
            if child.is_dir() and not child.is_symlink():
                shutil.rmtree(child)
            else:
                child.unlink()


for base in [
    FAMILY_ROOT / "examples" / "homo",
    FAMILY_ROOT / "validation" / "homo",
    ROOT / "tests",
]:
    remove_results(base)
    clear_out_dirs(base)
