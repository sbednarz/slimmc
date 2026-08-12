# Copolymer validation

Self-contained technical validation models used by `copo/Makefile`.

- `engine/`: seeded parser/engine validation models.
- `engine_channels/manifest.tsv`: black-box channel and pyslimmc-live groups.
- `run_channels.py`: checks or runs manifest entries with the shared family dispatcher.

The directory is part of the source package; Makefile targets must not point to files outside it.
