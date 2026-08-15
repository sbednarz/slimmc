# Releases and component tags

This repository is one monorepo with three independently versioned products:

| Component | Canonical version | Exact tag format | Release artifacts |
|---|---|---|---|
| Slimmc CLI + homo/copo engines | `VERSION` | `slimmc-v<VERSION>` | three native binary ZIP files |
| pyslimmc | `pyslimmc/_version.py` | `pyslimmc-v<VERSION>` | Python wheel and source distribution |
| pyslimmc-opt | `pyslimmc_opt/__init__.py` | `pyslimmc-opt-v<VERSION>` | Python wheel and source distribution |

The text after the component prefix and `v` must be byte-for-byte identical to
the version stored in the component files. Release workflows call
`scripts/check_release_tag.py` and stop before publication if the tag and files
differ. `make check-versions` also checks every build-time copy of all three
versions.

For the versions currently checked in, the only valid release tags are:

```text
slimmc-v5.0.2
pyslimmc-v5.0.0
pyslimmc-opt-v1.0.0
```

Python package tags must use the exact PEP 440 version spelling stored in the
component files; the checked-in spelling is the release contract.

## Publishing an existing checked-in version

First make sure the desired commit is on `main` and CI is green:

```bash
git switch main
git pull --ff-only
make check-versions
git status
```

Then create exactly one component tag. For example, to publish Slimmc:

```bash
git tag -a slimmc-v5.0.2 -m "Slimmc 5.0.2"
git push origin slimmc-v5.0.2
```

To publish the Python components from the same commit:

```bash
git tag -a pyslimmc-v5.0.0 -m "pyslimmc 5.0.0"
git push origin pyslimmc-v5.0.0

git tag -a pyslimmc-opt-v1.0.0 -m "pyslimmc-opt 1.0.0"
git push origin pyslimmc-opt-v1.0.0
```

Each tag creates a separate GitHub Release. A normal commit or push to `main`
does not publish a release.

## Changing a component version

Use the helper so all build-time copies change together:

```bash
python scripts/set_version.py slimmc 5.0.2
python scripts/set_version.py pyslimmc 5.0.0
python scripts/set_version.py pyslimmc-opt 1.0.0
```

Normally only the component being released is changed. Update its changelog,
README/release notes where relevant, then verify, commit and push:

```bash
make check-versions
git diff
git add -A
git commit -m "Prepare COMPONENT VERSION"
git push origin main
```

Wait for CI on `main` to pass. Only then create the matching component tag and
push it. Never reuse, move, or force-update a published release tag; increment
the release candidate or patch version instead.

## Slimmc native binaries

The `Release Slimmc binaries` workflow publishes:

| Archive suffix | Runtime contract | Intended systems |
|---|---|---|
| `linux-x86_64-glibc_2.28` | dynamically linked, glibc 2.28 or newer | current Linux and HPC clusters |
| `linux-x86_64-musl-static` | static ELF | portable Linux fallback |
| `windows-x86_64` | native 64-bit executable | Windows 10/11 x86-64 |

The glibc build uses a pinned Rocky Linux 8.10 builder image and requires the
builder's glibc to be exactly 2.28. The GitHub job itself remains on the Ubuntu
runner; no setup action is executed inside a container. Nim 2.2.10 is installed
inside the builder from the official Linux archive after checking its pinned
SHA-256. Both output executables are then inspected for imported GLIBC symbols.
The musl build is rejected if it contains an ELF interpreter or dynamic
`NEEDED` entries. All three variants run Slimmc version and homo/copo preflight
smoke tests before publication.

The workflows can also be run manually from GitHub Actions. Manual runs build
and retain artifacts but do not create a GitHub Release because they have no
component tag.

## See also

- [`DEVELOPMENT.md`](DEVELOPMENT.md) — Development workflow
- [`TESTING.md`](TESTING.md) — Testing and validation
- [`../reference/CLI.md`](../reference/CLI.md) — CLI reference
