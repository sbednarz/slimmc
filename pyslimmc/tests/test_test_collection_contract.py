from __future__ import annotations

import ast
from pathlib import Path


def test_every_test_module_declares_collectable_tests():
    root = Path(__file__).resolve().parent
    offenders: list[str] = []
    for path in sorted(root.glob("test_*.py")):
        if path.name == Path(__file__).name:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        collectable = any(
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_")
            or isinstance(node, ast.ClassDef) and node.name.startswith("Test")
            for node in tree.body
        )
        if not collectable:
            offenders.append(path.name)
    assert not offenders, "test modules with no pytest-collectable tests: " + ", ".join(offenders)
