# Homo versus effective one-monomer copo

This suite compares equivalent closed models executed by the standalone `homo`
and `copo` engines and read through public `pyslimmc`.

The copolymer grammar requires two or three declared monomers. Therefore the
"one-monomer copo" side is represented canonically by:

- reactive monomer `A`,
- spectator monomer `B` at exactly zero concentration,
- no channel that consumes or produces `B`.

The checker proves that `B` remains zero and that every copolymer chain satisfies
`dp == counts["A"]` and `counts["B"] == 0`.

Compared cases:

1. initiation;
2. initiation plus propagation;
3. build followed by isolated depropagation;
4. combination termination;
5. disproportionation termination.

Per-run stoichiometric identities are exact. Since the two engines have separate
implementations and may consume random numbers in a different order, kinetic
comparisons use ensemble means over several identical seed values rather than
requiring identical trajectories.

Run from repository root:

```bash
make test-homo-copo-equivalence
```
