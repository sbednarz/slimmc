# Homo phase C validation

Phase C validates the model-control layer without changing the model language:

- C11: exact `at`/`every` chronology, no implicit t=0 snapshot, same-time
  snapshot coalescing, and actions at `t_end`;
- C12: `set_k`, `add_k`, `set_temp`, `add_temp`, parameter-state snapshots,
  switching a channel off, and clean failure for a negative rate;
- C13: `set_c`/`add_c`, signed increments, monomer bookkeeping, no kinetic
  parameter-state increment, and clean failure below zero;
- C14: initial conditions, AND conditions, independent `when` lines, cascades
  caused by scheduled actions, one-shot behavior, `when_check_events`, and stop.

Run with:

```bash
cd copo
make test-phase-c
```
