from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
forbidden_dirs = {"__pycache__", "nimcache", "results", "opt_results", ".pytest_cache", "build", "dist"}
forbidden_suffixes = {".pyc", ".pyo"}
# Toolchains and virtual environments are dependencies, not build products of
# Slimmc.  In particular, setup-nim-action installs Nim below .nim_runtime in
# the GitHub Actions checkout; Nim itself contains a vendor/results directory.
# Skip these roots entirely while continuing to reject project results/ trees.
ignored_top_level_dirs = {".git", ".nim_runtime", ".venv", "venv"}
forbidden_binaries = {
    Path("common/tests/test_run_id"),
    Path("common/tests/test_run_id.exe"),
    Path("common/tests/test_model_contract"),
    Path("common/tests/test_model_contract.exe"),
    Path("common/tests/test_results_writers"),
    Path("common/tests/test_results_writers.exe"),
    Path("common/tests/test_storage_manifest"),
    Path("common/tests/test_storage_manifest.exe"),
    Path("homo/slimmc-resolved-test"),
    Path("homo/slimmc-resolved-test.exe"),
    Path("homo/slimmc-equivalence"),
    Path("homo/slimmc-equivalence.exe"),
}

bad = []
for path in ROOT.rglob("*"):
    relative = path.relative_to(ROOT)
    if relative.parts and relative.parts[0] in ignored_top_level_dirs:
        continue
    if path.is_dir() and path.name in forbidden_dirs:
        bad.append(path)
    # Committed example PNGs are documentation assets, not generated debris.
    elif path.is_dir() and path.name == "plots" and "examples" not in path.parts:
        bad.append(path)
    elif path.is_dir() and path.name.endswith(".egg-info"):
        bad.append(path)
    elif path.is_file() and path.suffix in forbidden_suffixes:
        bad.append(path)
    elif path.is_file() and path.parent == ROOT / "bin":
        bad.append(path)
    elif path.is_file() and relative in forbidden_binaries:
        bad.append(path)

if bad:
    raise SystemExit("Generated files remain:\n" + "\n".join(str(p.relative_to(ROOT)) for p in bad))
print("clean tree: PASS")
