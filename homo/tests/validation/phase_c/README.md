# Homo phase C validation

Phase C validates the model-control layer without changing the model language:

- H11: exact `at`/`every` chronology, periodic schedules include t=0, same-time
  snapshot coalescing, and actions at `t_end`;
- H12: `set_k`, `add_k`, `set_temp`, `add_temp`, parameter-state snapshots,
  switching a channel off, and clean failure for a negative rate;
- H13: `set_c`/`add_c`, signed increments, monomer bookkeeping, no kinetic
  parameter-state increment, and clean failure below zero;
- H14: initial conditions, AND conditions, independent `when` lines, cascades
  caused by scheduled actions, one-shot behavior, `when_check_events`, and stop.

Run with:

```bash
cd homo
make test-phase-c
```
