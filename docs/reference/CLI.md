# CLI

The family dispatcher always uses this order:

```text
slimmc [options] model.model
```

```bash
slimmc model.model
slimmc --check model.model
slimmc --debug model.model
slimmc --trace-channels 10000 model.model
slimmc --output-root results model.model
```

The model file is always the final argument. The dispatcher recognizes homo
and copo models from model contents.

Information commands:

```bash
slimmc             # component versions and short usage
slimmc -h          # short model-running guide
slimmc --help      # same as -h
slimmc --version   # component versions and available build information
slimmc -v          # same as --version
slimmc -V          # same as --version
```

The version report shows the single Slimmc version shared by the family CLI,
homo engine, and copo engine, followed by the independent pyslimmc and
pyslimmc-opt versions. Build information includes the build mode, optimization,
Nim version, compiler backend, target platform, compilation time, and optional
Git provenance when it was embedded during compilation.

## Run summaries

Native binary releases also contain `slimmc-summary`, a read-only command for a
Slimmc Storage run directory:

```bash
slimmc-summary RUN
slimmc-summary RUN --format json
slimmc-summary RUN -o summary.txt
slimmc-summary RUN --format json -o summary.json
slimmc-summary --version
```

Text output includes run/engine identity, lifecycle and validation status, seed,
final time and KMC event, snapshot count, final DP/molar-mass moments, peak
memory and result-directory size. JSON mode exposes the same compact summary for
scripts. The command reads an existing run; it does not modify Storage.

## Common model/check failures

Start with `slimmc --check model.model`; parser/validator diagnostics include
model-line context where available. Common causes are:

- an unknown monomer/species/pool name in a reaction or condition;
- a rate name used without a matching `rate` declaration;
- a non-positive or otherwise invalid `kmc_volume`, time, rate, or action value;
- a homo `when X ...` condition that does not use the required explicit monomer
  form;
- an `output_dir` path that violates the model-path contract (for example `..`,
  spaces, glob characters, or invalid path segments);
- engine-specific syntax copied from homo into copo or vice versa;
- requesting chain-derived analysis from a snapshot saved with `save` but not
  `save_chains`; this is a data-availability issue rather than a parser error.

A warning from `--check` is not automatically fatal. In particular,
discretization warnings should trigger a review of `V_kMC` and convergence, not
a mechanical increase until the warning disappears. Exact syntax is in
[`HOMO.md`](HOMO.md), [`COPO.md`](COPO.md), and the shared
[`../MODEL_SYNTAX.md`](../MODEL_SYNTAX.md).

## See also

- [`HOMO.md`](HOMO.md) — Homo engine reference
- [`COPO.md`](COPO.md) — Copo engine reference
- [`../QUICKSTART.md`](../QUICKSTART.md) — Quick start
