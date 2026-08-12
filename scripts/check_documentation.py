#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
CANONICAL = [ROOT / "README.md", ROOT / "AGENTS.md", *sorted(DOCS.rglob("*.md"))]
FORBIDDEN = (
    "run.history", "pyslimmc.homo", "pyslimmc.copo", "dlogM",
    "legacy_dlog", "legacy gaussian", "basis: number | mass | z",
    "run.runinfo", "run.parameter_states", "run.oligomers",
)


def markdown_links(path: Path) -> list[str]:
    return re.findall(r"(?<!!)\[[^\]]+\]\(([^)#]+)(?:#[^)]+)?\)", path.read_text(encoding="utf-8"))


def check_links() -> list[str]:
    errors: list[str] = []
    for path in CANONICAL:
        for target in markdown_links(path):
            if "://" in target or target.startswith("mailto:"):
                continue
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                errors.append(f"broken link: {path.relative_to(ROOT)} -> {target}")
    return errors


def check_forbidden() -> list[str]:
    errors: list[str] = []
    for path in CANONICAL:
        text = path.read_text(encoding="utf-8")
        for token in FORBIDDEN:
            if token in text:
                errors.append(f"forbidden API name in {path.relative_to(ROOT)}: {token}")
    return errors


def check_versions() -> list[str]:
    errors: list[str] = []
    sys.path.insert(0, str(ROOT))
    try:
        expected = importlib.import_module("pyslimmc").__version__
    finally:
        sys.path.pop(0)
    for path in (ROOT / "README.md", DOCS / "reference" / "PYSLIMMC_API.md"):
        if expected not in path.read_text(encoding="utf-8"):
            errors.append(f"missing pyslimmc {expected} in {path.relative_to(ROOT)}")
    return errors


def check_public_exports() -> list[str]:
    sys.path.insert(0, str(ROOT))
    try:
        packages = (
            (importlib.import_module("pyslimmc"), DOCS / "reference" / "PYSLIMMC_API.md"),
            (importlib.import_module("pyslimmc_opt"), DOCS / "reference" / "PYSLIMMC_OPT.md"),
        )
    finally:
        sys.path.pop(0)
    errors: list[str] = []
    for pkg, path in packages:
        reference = path.read_text(encoding="utf-8")
        errors.extend(
            f"public export missing from {path.relative_to(ROOT)}: {pkg.__name__}.{name}"
            for name in pkg.__all__
            if name not in reference
        )
    return errors


def _of_tokens(text: str, start: str, end: str, indent: int) -> set[str]:
    section = text.split(start, 1)[1].split(end, 1)[0]
    pattern = "^" + (" " * indent) + r'of\s+((?:"[^"]+"(?:,\s*)?)+):'
    found: set[str] = set()
    for group in re.findall(pattern, section, flags=re.M):
        found.update(re.findall(r'"([^"]+)"', group))
    return found


def check_language_reference() -> list[str]:
    errors: list[str] = []
    specifications = (
        (
            ROOT / "homo" / "src" / "slimmc_parser.nim",
            DOCS / "reference" / "HOMO.md",
            (
                ("proc parseParam", "proc parseSpeciesDecl", 2),
                ("proc parseActionKind", "proc parseDesc", 2),
                ("proc parseMacroLine", "proc parseReactionLine", 2),
                ("proc parseModel*(path", "if m.V <= 0.0", 4),
            ),
        ),
        (
            ROOT / "copo" / "src" / "copo_parser.nim",
            DOCS / "reference" / "COPO.md",
            (
                ("proc parseActionKind", "proc validateActionArgs", 2),
                ("proc parseMacro", "proc inferPoolMetadata", 2),
                ("case p", '    of "desc"', 6),
                ("case key", "require(result.V > 0", 4),
            ),
        ),
    )
    ignored = {
        "volume",       # rejected migration diagnostic, not valid syntax
        "memory_limit", # intentionally not part of the public documentation; use at_memory
    }
    for source_path, doc_path, sections in specifications:
        source = source_path.read_text(encoding="utf-8")
        reference = doc_path.read_text(encoding="utf-8")
        tokens: set[str] = set()
        for start, end, indent in sections:
            tokens.update(_of_tokens(source, start, end, indent))
        for token in sorted(tokens - ignored):
            if not re.search(
                rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])",
                reference,
            ):
                errors.append(
                    f"parser token missing from {doc_path.name}: {token} "
                    f"(source {source_path.relative_to(ROOT)})"
                )
    defaults = ("298.15", "10000000000", "12345", "2147483647", "repeat_units", "composition")
    for doc_path in (DOCS / "reference" / "HOMO.md", DOCS / "reference" / "COPO.md"):
        reference = doc_path.read_text(encoding="utf-8")
        for value in defaults:
            if value not in reference:
                errors.append(f"shared parser default missing from {doc_path.name}: {value}")
    return errors


def check_signature_inventory() -> list[str]:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "update_api_signatures.py"), "--check"],
        text=True,
        capture_output=True,
    )
    return [] if proc.returncode == 0 else [proc.stdout.strip() or proc.stderr.strip()]


def check_integration_coverage() -> list[str]:
    reference = (DOCS / "development" / "INTEGRATION_COVERAGE.md").read_text(encoding="utf-8")
    models = sorted((ROOT / "tests" / "integration" / "models").glob("*.model"))
    return [
        f"integration model missing from coverage matrix: {model.name}"
        for model in models
        if f"`{model.name}`" not in reference
    ]


def check_generated_toc() -> list[str]:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "update_markdown_toc.py"),
         "--check", str(DOCS / "development" / "DEVELOPMENT.md")],
        text=True, capture_output=True,
    )
    return [] if proc.returncode == 0 else [proc.stdout.strip() or proc.stderr.strip()]


def find_engine(cli_arg: str | None) -> Path | None:
    candidates = [cli_arg, os.environ.get("SLIMMC_CLI"), str(ROOT / "bin" / "slimmc")]
    for raw in candidates:
        if raw and Path(raw).is_file():
            return Path(raw).resolve()
    return None


def check_model_blocks(engine: Path | None, require_engine: bool) -> list[str]:
    if engine is None:
        if require_engine:
            return ["model checks requested but no CLI was provided; use --engine or SLIMMC_CLI"]
        print("documentation model checks: SKIP (provide --engine or SLIMMC_CLI)")
        return []
    errors: list[str] = []
    n = 0
    for doc in (DOCS / "QUICKSTART.md", DOCS / "COOKBOOK.md", DOCS / "reference" / "HOMO.md", DOCS / "reference" / "COPO.md"):
        text = doc.read_text(encoding="utf-8")
        blocks = re.findall(r"```(?:text|model)\n(.*?)\n```", text, flags=re.S)
        for block in blocks:
            if "param output_dir" not in block or "monomer " not in block:
                continue
            n += 1
            with tempfile.TemporaryDirectory() as td:
                model = Path(td) / f"doc_{n}.model"
                model.write_text(block + "\n", encoding="utf-8")
                proc = subprocess.run([str(engine), "--check", str(model)], text=True, capture_output=True)
                if proc.returncode:
                    errors.append(f"model block failed in {doc.name}:\n{proc.stdout}{proc.stderr}")
    print(f"documentation model checks: {n} block(s) PASS" if not errors else f"documentation model checks: {n} checked")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", help="path to compiled unified slimmc CLI")
    parser.add_argument("--require-engine", action="store_true", help="fail instead of skipping model checks")
    args = parser.parse_args()

    errors = []
    errors += check_links()
    errors += check_forbidden()
    errors += check_versions()
    errors += check_public_exports()
    errors += check_language_reference()
    errors += check_signature_inventory()
    errors += check_integration_coverage()
    errors += check_generated_toc()
    errors += check_model_blocks(find_engine(args.engine), args.require_engine)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("documentation contract: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
