#!/usr/bin/env python3
"""Guard the release workflow against host/container path regressions."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(text: str, needle: str, source: Path) -> None:
    if needle not in text:
        raise SystemExit(f"{source}: missing required release contract: {needle!r}")


def reject(text: str, needle: str, source: Path) -> None:
    if needle in text:
        raise SystemExit(f"{source}: forbidden release configuration: {needle!r}")


def main() -> None:
    workflow_path = ROOT / ".github/workflows/release.yml"
    dockerfile_path = ROOT / ".github/release/linux-glibc-2.28.Dockerfile"
    build_script_path = ROOT / "scripts/ci/build_linux_glibc_2_28.sh"

    workflow = workflow_path.read_text(encoding="utf-8")
    start = workflow.index("  linux-glibc:\n")
    end = workflow.index("  linux-musl:\n", start)
    job = workflow[start:end]

    reject(job, "container:", workflow_path)
    reject(job, "setup-nim-action", workflow_path)
    reject(job, "setup-python", workflow_path)
    require(job, "docker build --pull", workflow_path)
    require(job, ".github/release/linux-glibc-2.28.Dockerfile", workflow_path)
    require(job, "scripts/ci/build_linux_glibc_2_28.sh", workflow_path)
    require(job, "if-no-files-found: error", workflow_path)

    dockerfile = dockerfile_path.read_text(encoding="utf-8")
    require(dockerfile, "FROM rockylinux/rockylinux:8.10", dockerfile_path)
    require(dockerfile, "ARG NIM_VERSION=2.2.10", dockerfile_path)
    require(
        dockerfile,
        "0a3a38752e97e9d44aa479b3a7b37336dfe0176daf22ee5b5218ad0991ecd211",
        dockerfile_path,
    )
    require(dockerfile, "sha256sum -c -", dockerfile_path)
    require(dockerfile, 'GNU_LIBC_VERSION', dockerfile_path)
    require(dockerfile, "findutils", dockerfile_path)
    require(dockerfile, "python39", dockerfile_path)

    build_script = build_script_path.read_text(encoding="utf-8")
    require(build_script, 'readonly EXPECTED_GLIBC="2.28"', build_script_path)
    require(build_script, 'uname -m', build_script_path)
    require(build_script, 'actual_glibc" != "$EXPECTED_GLIBC', build_script_path)
    require(build_script, "--glibc-max", build_script_path)
    require(build_script, "FRP_ALLCHAN01.model", build_script_path)
    require(build_script, "C01_init.model", build_script_path)
    require(build_script, "check_build_provenance.py", build_script_path)
    require(build_script, "--require-git", build_script_path)
    require(build_script, "package_binary_release.py", build_script_path)

    print("release configuration: PASS (isolated Rocky Linux 8.10/glibc 2.28 build)")


if __name__ == "__main__":
    main()
