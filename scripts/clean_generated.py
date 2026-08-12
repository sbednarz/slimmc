#!/usr/bin/env python3
"""Remove family-level outputs that are reproducible from source files."""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def remove_tree(path: Path) -> None:
    if path.is_dir():
        print("removing", path.relative_to(ROOT))
        shutil.rmtree(path)


for path in (
    ROOT / "output", ROOT / "tmp", ROOT / ".pytest_cache", ROOT / "artifacts",
    ROOT / "build", ROOT / "dist", ROOT / "pyslimmc_opt" / "build",
    ROOT / "pyslimmc_opt" / "dist",
):
    remove_tree(path)

literature = ROOT / "literature"
if literature.exists():
    for path in sorted(literature.rglob("results*"), reverse=True):
        remove_tree(path)
    for out_dir in sorted(literature.rglob("out"), reverse=True):
        if not out_dir.is_dir():
            continue
        for child in list(out_dir.iterdir()):
            if child.name in {".keep", ".gitkeep"}:
                continue
            if child.is_dir() and not child.is_symlink():
                remove_tree(child)
            else:
                print("removing", child.relative_to(ROOT))
                child.unlink()

simulations = ROOT / "simulations"
if simulations.exists():
    for path in sorted(simulations.glob("*/results*")):
        remove_tree(path)

examples = ROOT / "examples"
if examples.exists():
    # plots/ contains committed reference PNGs embedded in example READMEs.
    # Only reproducible simulation and optimization storage is removed.
    for name in ("results", "opt_results"):
        for path in sorted(examples.rglob(name), reverse=True):
            remove_tree(path)

# Black-box validation writes run-local Storage trees beside its model files.
for component in ("homo", "copo"):
    tests = ROOT / component / "tests"
    if tests.exists():
        for path in sorted(tests.rglob("results"), reverse=True):
            remove_tree(path)

for path in sorted(ROOT.rglob("*.egg-info"), reverse=True):
    remove_tree(path)

# Component Makefiles clean their own caches. This family-level pass also
# catches caches produced by shared tests, examples, validation and reports.
for path in sorted(ROOT.rglob("__pycache__"), reverse=True):
    remove_tree(path)
for path in sorted(ROOT.rglob(".pytest_cache"), reverse=True):
    remove_tree(path)
for path in sorted(ROOT.rglob(".ipynb_checkpoints"), reverse=True):
    remove_tree(path)
for path in ROOT.rglob("*.pyc"):
    if path.is_file():
        print("removing", path.relative_to(ROOT))
        path.unlink()
