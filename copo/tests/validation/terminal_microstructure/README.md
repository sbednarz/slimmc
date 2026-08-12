# Copo terminal, penultimate and microstructure validation

This group complements the low-level chemistry tests with full black-box runs.
It verifies all terminal and penultimate propagation channels, pool metadata,
composition, engine dyads/triads/block counts, reconstruction from literal
sequences in `full` mode, and parity of chemistry and aggregate microstructure
between `sequence_mode full` and `sequence_mode composition`.

Run with:

    make test-terminal-microstructure
