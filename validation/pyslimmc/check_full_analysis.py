import numpy as np

from pyslimmc.full_analysis import sequence_stats, block_lengths, transition_matrix, microstructure_by_dp

class Composition:
    names = ("IA", "BMA")

class FakeChains:
    def __init__(self):
        self.sequences = (
            ("IA", "IA", "BMA", "BMA", "IA"),
            ("BMA", "IA", "BMA"),
            ("IA", "IA", "IA"),
        )
        self.count = np.array([2, 3, 5], dtype=np.uint64)
        self.dp = np.array([5, 3, 3], dtype=np.uint64)
        self.composition = Composition()
        self.has_sequences = True

c = FakeChains()
s = sequence_stats(c, progress=False)
assert s.transition_count.tolist() == [2, 2, 0]
assert np.allclose(s.transition_fraction, [0.5, 1.0, 0.0])
assert s.block_count["IA"].tolist() == [2, 1, 1]
assert s.max_block_length["IA"].tolist() == [2, 1, 3]
assert np.allclose(s.mean_block_length["IA"], [1.5, 1.0, 3.0])

b = block_lengths(c, "IA", progress=False)
# IA blocks: row1 lengths 2 and 1, each weighted 2; row2 length 1 weighted 3; row3 length 3 weighted 5
assert b.length.tolist() == [1, 2, 3]
assert np.allclose(b.count, [5, 2, 5])
assert np.isclose(b.fraction.sum(), 1.0)

m = transition_matrix(c, normalize=None, progress=False)
# weighted transitions: IA->IA 2 + 10, IA->BMA 2 + 3, BMA->IA 2 + 3, BMA->BMA 2
assert np.allclose(m.values, [[12, 5], [5, 2]])
mr = transition_matrix(c, normalize="row", progress=False)
assert np.allclose(mr.values.sum(axis=1), [1, 1])

r = microstructure_by_dp(c, "transition_fraction", bins=[2.5, 3.5, 5.5], progress=False)
assert r.chain_count.tolist() == [8.0, 2.0]
assert np.allclose(r.mean, [3/8, 0.5])
print("full analysis: PASS")
