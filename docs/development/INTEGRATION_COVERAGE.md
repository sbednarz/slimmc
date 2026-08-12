# Integration coverage

The independent suite contains 92 tests in `tests/integration/`. It executes
compact private fixtures through the unified CLI and reads only real Slimmc
Storage through public `pyslimmc`/`pyslimmc-opt` APIs. It never reads or runs

## Engine scenarios

| Scenario | Engine | Contract covered |
|---|---|---|
| `homo.model` | homo | initiation, propagation, disproportionation, feed, temperature action, snapshots, chains, moments, MWD/CLD, balances |
| `copo.model` | copo | terminal binary propagation, combination, feed, full sequences, composition, microstructure, spectra |
| `homo_mechanisms.model` | homo | `rxn`, every homo `macro` kind, Arrhenius, all scheduled mutable/output actions, all `var` target kinds, `with_end_groups`, memory declarations |
| `copo_mechanisms.model` | copo | `rxn`, every copo `macro` kind, Arrhenius, all scheduled mutable/output actions, feed, `with_end_groups`, exact sequences |
| `homo_stop.model` | homo | multi-line `when`, explicit save/save_chains, clean conditional `stop` |
| `copo_penultimate.model` | copo | 12-channel explicit binary penultimate topology and penultimate analysis |
| `copo_terpolymer.model` | copo | three monomers, nine propagation transitions, ternary compositions and binary-analysis rejection |
| `copo_composition.model` | copo | `sequence_mode composition`, aggregate dyads, sequence-only analysis rejection |
| `homo_failed.model` | homo | genuine runtime failure and diagnostic incomplete Storage |
| `homo_interrupted.model` | homo | genuine SIGINT handling and interrupted Storage |
| `opt.model` | homo/opt | three-trial `t_end` optimization, resume, artifacts and best-model rerun |
| `opt_surface.model` | homo/opt | fixed-rate and monomer substitution plus custom `build` and independent verification |

`at_memory` crosses the parser/dispatcher/Storage pipeline in
the integration models. Deliberately reaching a machine-memory threshold is
kept in targeted technical tests because it is environment-dependent.

## Public behavior matrix

| Public area | Integration assertions |
|---|---|
| CLI | help/version aliases, `--check`, dispatcher selection, missing/invalid model failures |
| Lifecycle | completed, failed, interrupted, strict direct open, diagnostic opt-in, `Runs` status namespaces |
| Chronology | snapshot IDs, time/event axes, no artificial `t=0`, final/last, nearest selectors |
| State/process | counts, concentrations, physical/KMC volume, initial amounts, feed events, cumulative doses, balances |
| Kinetics/actions | channel firings, productive intervals, Arrhenius evolution, `set_k`, `add_k`, `set_temp`, `add_temp`, `set_c`, `add_c` |
| Chains/masses | live/dead counts, compressed records, sequences, pools/origins, repeat-unit and end-group mass audits |
| Distributions | moments, MWD, CLD, neutral chain-mass spectrum, chain counts |
| Copolymer analysis | `f`, `F.ins`, `F.int`, `F.cum`, terminal, penultimate, terpolymer and microstructure applicability |
| Collections | scan, deterministic selection, status filters and lifecycle namespaces |
| Optimization | real trials, seeds, GP/EI budget, fixed-rate/monomer substitution, custom `build`, failure, resume guard, best-model rerun and independent verification |

Parser rejection rules, stochastic laws, numerical boundaries, Storage schema
invariants, and chemistry-specific regression cases remain in unit,
black-box, and `validation/` suites. Integration tests establish that the
public layers work together; they do not replace those more focused tests.

## See also

- [`TESTING.md`](TESTING.md) — Testing strategy
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — Architecture
- [`../reference/STORAGE.md`](../reference/STORAGE.md) — Storage format
