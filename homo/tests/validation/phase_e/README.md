# Homo validation phase E

Phase E validates user-visible results rather than new reaction chemistry:

- H20: exact chain molar masses for `repeat_units` and `with_end_groups`;
- H21: deterministic aggregation of identical chain records and contiguous IDs;
- H22: state-only versus chain-bearing snapshots, no synthetic time zero, and coalescing of `save` with `save_chains`;
- H23: completed, failed, interrupted, zero-propensity and explicit-stop finalization contracts.

The checker uses the public `pyslimmc` API for semantic checks and reads canonical storage arrays only for low-level format invariants that are not part of the public analysis API.
