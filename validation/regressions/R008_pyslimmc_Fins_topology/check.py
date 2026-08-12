from pathlib import Path
import numpy as np
import pyslimmc as sl
run=sl.open(Path(__file__).parent/'results/main')
fa=np.asarray(run.F.ins['A'],dtype=float)
assert np.all(np.isfinite(fa))
assert np.any(np.abs(fa-np.asarray(run.f['A'],dtype=float))>1e-4)
print('R008 PASS', fa[-1])
