from pathlib import Path
import pyslimmc as sl
run=sl.open(Path(__file__).parents[1]/'R008_pyslimmc_Fins_topology/results/main')
assert run.snapshots_with_chains
assert run.first_with_chains.has_chains
assert run.last_with_chains.has_chains
chains=run.last_with_chains.chains
assert chains.n_records == len(chains)
assert chains.n_chains == int(chains.count.sum())
print('R009 PASS', len(run.snapshots_with_chains), chains.n_records, chains.n_chains)
