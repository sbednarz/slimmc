# Copo phase E validation

Black-box validation through the shared `slimmc` dispatcher and public
`pyslimmc` API.

- C20: repeat-unit and end-group molar masses for binary copolymers.
- C21: aggregation keys retain composition, terminal metadata and sequences.
- C22: state-only versus chain snapshots, `F.cum` for every snapshot, and
  coalescing of `save` plus `save_chains` at the same time.
- C23: clean completion, zero-propensity completion, conditional stop,
  controlled failure and SIGINT interruption.

Run from `copo/`:

    make test-phase-e
