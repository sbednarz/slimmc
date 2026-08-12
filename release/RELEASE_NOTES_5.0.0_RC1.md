# Slimmc 5.0.0-rc.1

Release-candidate versions:

- Slimmc 5.0.0-rc.1 — one version shared by the CLI, homo engine, and copo
  engine;
- pyslimmc 4.0.0rc1;
- pyslimmc-opt 1.0.0rc1.

The public command-line contract is:

```text
slimmc [options] model.model
```

Bare `slimmc` prints the three public component versions and short usage.
`slimmc --version` additionally prints available compilation details.
`slimmc -h` and `slimmc --help` show a short model-running guide.

The model family continues to be selected internally from the model contents,
but this implementation detail is not printed by the information commands.
Slimmc Storage remains at format version 1.2.0.
