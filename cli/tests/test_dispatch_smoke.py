#!/usr/bin/env python3
"""End-to-end smoke test for cli/slimmc_cli.nim (the slimmc dispatcher),
run against real model files already present in this package.

This is the "Tier 2" end-to-end dispatcher check referenced in
docs/development/TESTING.md ("End-to-end checks additionally inspect the resulting
engine identity"). It complements, and does not replace,
cli/tests/test_dispatch.nim's pure dispatch-logic unit tests (which test
dispatch_logic.nim in isolation, without linking either engine).

Usage: test_dispatch_smoke.py <path-to-slimmc-binary>
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def run(binary: Path, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([str(binary), *args], cwd=cwd, capture_output=True, text=True)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_dispatch_smoke.py <path-to-slimmc-binary>")
    binary = Path(sys.argv[1]).resolve()
    require(binary.exists(), f"{binary} not found -- run `make build` first")

    homo_model = ROOT / "homo" / "tests" / "regression" / "frp_all_channels_seeded" / "FRP_ALLCHAN01.model"
    copo_model = ROOT / "examples" / "basic" / "064_ideal_binary_copo" / "model.model"
    require(homo_model.exists(), f"fixture model missing: {homo_model}")
    require(copo_model.exists(), f"fixture model missing: {copo_model}")

    # --- no arguments prints all component versions and short usage ---
    result = run(binary)
    require(result.returncode == 0, result.stderr)
    slimmc_version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    pyslimmc_version = (ROOT / "pyslimmc" / "_version.py").read_text(encoding="utf-8").split('"')[1]
    opt_version = (ROOT / "pyslimmc_opt" / "__init__.py").read_text(encoding="utf-8").split('"')[-2]
    for expected in (
        f"Slimmc {slimmc_version}",
        f"pyslimmc {pyslimmc_version}",
        f"pyslimmc-opt {opt_version}",
    ):
        require(expected in result.stdout, f"missing version line {expected!r}:\n{result.stdout}")
    require("Usage: slimmc [options] model.model" in result.stdout, result.stdout)
    print("[OK] bare slimmc prints all component versions and short usage")

    # --- --version adds compilation information ---
    result = run(binary, "--version")
    require(result.returncode == 0, result.stderr)
    for expected in ("Build", "mode:", "optimization:", "Nim:", "backend:", "target:", "compiled:"):
        require(expected in result.stdout, f"missing build field {expected!r}:\n{result.stdout}")
    print("[OK] slimmc --version prints component and build information")

    # --- -h/--help show the new options-before-model contract ---
    for help_flag in ("-h", "--help"):
        result = run(binary, help_flag)
        require(result.returncode == 0, result.stderr)
        require("slimmc [options] model.model" in result.stdout, result.stdout)
        require("--check model.model" in result.stdout, result.stdout)
    print("[OK] -h and --help show a short model-running guide")

    # --- routes a 1-monomer model to homo, in-process ---
    result = run(binary, "--check", str(homo_model))
    require(result.returncode == 0, f"homo --check failed: {result.stderr}")
    require("scheduled actions" in result.stdout, "expected homo's --check output, got:\n" + result.stdout)
    print("[OK] slimmc routes a 1-monomer model to homo and runs --check in-process")

    # --- routes a 2/3-monomer model to copo, in-process ---
    result = run(binary, "--check", str(copo_model))
    require(result.returncode == 0, f"copo --check failed: {result.stderr}")
    require("discretization preflight" in result.stdout.lower(),
            "expected copo's --check output, got:\n" + result.stdout)
    print("[OK] slimmc routes a 2-3 monomer model to copo and runs --check in-process")

    # --- 0 monomers is a legitimate homo case (pure-kinetics model, no
    #     polymer chain growth at all). Use a real
    #     model from homo's own validation suite (20 of its 44 engine
    #     validation cases have zero monomer declarations), not a
    #     synthetic one, since a synthetic minimal file would also be
    #     missing other required fields and fail for unrelated reasons. ---
    zero_monomer_model = ROOT / "validation" / "regressions" / "R001_homo_termination_tend" / "model.model"
    require(zero_monomer_model.exists(), f"fixture model missing: {zero_monomer_model}")
    result = run(binary, "--check", str(zero_monomer_model))
    require(result.returncode == 0, f"0-monomer pure-kinetics model should route to homo: {result.stderr}")
    require("scheduled actions" in result.stdout, "expected homo's --check output, got:\n" + result.stdout)
    print("[OK] slimmc routes a 0-monomer pure-kinetics model to homo (not a dispatcher error)")

    # --- dispatcher-level rejection still happens before either engine
    #     runs, for a genuinely unroutable monomer count (>3) ---
    with tempfile.TemporaryDirectory() as tmp:
        bad_model = Path(tmp) / "too_many_monomers.model"
        bad_model.write_text(
            "monomer A 0.2 100.0\nmonomer B 0.2 100.0\n"
            "monomer C 0.2 100.0\nmonomer D 0.2 100.0\n"
        )
        result = run(binary, str(bad_model))
        require(result.returncode != 0, "a model with >3 monomers must be rejected")
        require("monomer` declarations found" in result.stderr, result.stderr)
    print("[OK] slimmc rejects a >3-monomer model at the dispatcher level")

    # --- an engine-level error (bad flag) surfaces with the real engine's
    #     own usage/error text, forwarded through unchanged ---
    result = run(binary, "--this-flag-does-not-exist", str(copo_model))
    require(result.returncode != 0, "an unknown CLI flag must fail")
    require("unknown option: --this-flag-does-not-exist" in result.stderr, result.stderr)
    print("[OK] slimmc rejects an unknown CLI option")

    # --- old options-after-model syntax is rejected clearly ---
    result = run(binary, str(homo_model), "--check")
    require(result.returncode != 0, "old model-before-options syntax must fail")
    require("model file must be last" in result.stderr, result.stderr)
    print("[OK] slimmc enforces options before the final model file")

    # --- full (non-check) run actually completes end to end ---
    with tempfile.TemporaryDirectory() as tmp:
        local_model = Path(tmp) / homo_model.name
        local_model.write_text(homo_model.read_text())
        result = run(binary, str(local_model), cwd=Path(tmp))
        require(result.returncode == 0, f"full homo run failed: {result.stderr}\n{result.stdout}")
        require((Path(tmp) / "results").is_dir(), "expected a results/ directory from a full run")
    print("[OK] slimmc completes a full (non-check) simulation end to end")

    print("[test_dispatch_smoke] ok")


if __name__ == "__main__":
    main()
