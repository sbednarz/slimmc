from pathlib import Path
import pyslimmc as sl
run=sl.open(Path(__file__).parents[1]/'R006_copo_terpoly_feed_species_padding/results/main')
feed=next(iter(run.feeds.values()))
assert len(feed.events)==feed.events.n_events==feed.events.time.size
print('R009 feed PASS', len(feed.events))
