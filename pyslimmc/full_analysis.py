from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .core import DataUnavailableError
from .options import options
from .composition_analysis import _axes, _save, _bin_edges, _readonly


def _progress(name: str, enabled: bool, done: int, total: int, state: dict) -> None:
    if not enabled or total <= 0:
        return
    percent = int(done * 100 / total)
    for mark in (25, 50, 75):
        if percent >= mark and mark not in state:
            print(f"[{name}] {mark}%")
            state[mark] = True


@dataclass(frozen=True)
class SequenceStats:
    names: tuple[str, ...]
    transition_count: np.ndarray
    transition_fraction: np.ndarray
    block_count: Mapping[str, np.ndarray]
    max_block_length: Mapping[str, np.ndarray]
    mean_block_length: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class BlockLengthDistribution:
    monomer: str | None
    length: np.ndarray
    count: np.ndarray
    fraction: np.ndarray

    def plot(self, path=None, *, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        ax.plot(self.length, self.count, marker="o")
        ax.set_xlabel("block length")
        ax.set_ylabel("block count")
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class TransitionMatrix:
    names: tuple[str, ...]
    values: np.ndarray
    normalization: str | None

    def plot(self, path=None, *, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        image = ax.imshow(self.values, origin="upper", aspect="equal")
        ax.set_xticks(np.arange(len(self.names)), self.names)
        ax.set_yticks(np.arange(len(self.names)), self.names)
        ax.set_xlabel("next monomer"); ax.set_ylabel("current monomer")
        fig.colorbar(image, ax=ax, label="transition fraction" if self.normalization else "transition count")
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class MotifAxisMap:
    axis_name: str
    edges: np.ndarray
    motifs: tuple[str, ...]
    count: np.ndarray
    fraction: np.ndarray

    @property
    def center(self) -> np.ndarray:
        return (self.edges[:-1] + self.edges[1:]) / 2.0

    def plot(self, path=None, *, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        image = ax.imshow(self.fraction.T, origin="lower", aspect="auto",
                          extent=(self.edges[0], self.edges[-1], -0.5, len(self.motifs)-0.5),
                          vmin=0, vmax=np.nanmax(self.fraction))
        ax.set_yticks(np.arange(len(self.motifs)), self.motifs)
        ax.set_xlabel(self.axis_name); ax.set_ylabel("ordered motif")
        fig.colorbar(image, ax=ax, label="motif fraction")
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


def _motif_axis_map(chains, axis, edges, order: int, axis_name: str) -> MotifAxisMap:
    names = tuple(chains.composition.names)
    motifs = tuple("".join(parts) for parts in __import__("itertools").product(names, repeat=order))
    index = {motif: i for i, motif in enumerate(motifs)}
    counts = np.zeros((len(edges)-1, len(motifs)), dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    sequences = chains.sequences
    for row, (seq, weight, value) in enumerate(zip(sequences, weights, axis)):
        b = int(np.searchsorted(edges, value, side="right") - 1)
        if value == edges[-1]: b = len(edges)-2
        if b < 0 or b >= len(edges)-1: continue
        for i in range(max(0, len(seq)-order+1)):
            motif = "".join(seq[i:i+order])
            counts[b, index[motif]] += weight
    totals = counts.sum(axis=1, keepdims=True)
    fractions = np.divide(counts, totals, out=np.full_like(counts, np.nan), where=totals != 0)
    return MotifAxisMap(axis_name, _readonly(edges), motifs, _readonly(counts), _readonly(fractions))


def dyads_by_dp(chains, *, bins=16) -> MotifAxisMap:
    dp = np.asarray(chains.dp, dtype=float)
    if dp.size == 0: raise DataUnavailableError("chain population is empty")
    edges = _bin_edges(dp, bins, integer=True)
    return _motif_axis_map(chains, dp, edges, 2, "DP")


def triads_by_composition(chains, monomer: str, *, bins=12) -> MotifAxisMap:
    if monomer not in chains.composition.names: raise KeyError(monomer)
    fraction = np.asarray(chains.composition.fractions[monomer], dtype=float)
    edges = np.linspace(0.0, 1.0, int(bins)+1) if isinstance(bins, int) else np.asarray(bins, dtype=float)
    return _motif_axis_map(chains, fraction, edges, 3, f"{monomer} fraction in chain")


@dataclass(frozen=True)
class MicrostructureByDP:
    statistic: str
    monomer: str | None
    dp_left: np.ndarray
    dp_right: np.ndarray
    dp_center: np.ndarray
    record_count: np.ndarray
    chain_count: np.ndarray
    mean: np.ndarray
    median: np.ndarray
    q25: np.ndarray
    q75: np.ndarray

    def plot(self, path=None, *, statistic="mean", interval=None, style="screen", ax=None,
             span=None, dpi=300, title=None):
        if statistic not in {"mean", "median"}:
            raise ValueError("statistic must be 'mean' or 'median'")
        if interval not in {None, "iqr"}:
            raise ValueError("interval must be None or 'iqr'")
        fig, ax = _axes(ax, style=style, span=span, title=title)
        values = getattr(self, statistic)
        ax.plot(self.dp_center, values)
        if interval == "iqr":
            ax.fill_between(self.dp_center, self.q25, self.q75, alpha=0.18)
        ax.set_xlabel("DP")
        label = self.statistic if self.monomer is None else f"{self.statistic} ({self.monomer})"
        ax.set_ylabel(label.replace("_", " "))
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


def _sequence_rows(chains):
    if not chains.has_sequences:
        raise DataUnavailableError("complete stored sequences are required")
    return chains.sequences


def _slice_sequence_stats(stats: SequenceStats, indices: np.ndarray) -> SequenceStats:
    indices = np.asarray(indices, dtype=np.int64)
    return SequenceStats(
        names=stats.names,
        transition_count=_readonly(stats.transition_count[indices]),
        transition_fraction=_readonly(stats.transition_fraction[indices]),
        block_count=MappingProxyType({k: _readonly(v[indices]) for k, v in stats.block_count.items()}),
        max_block_length=MappingProxyType({k: _readonly(v[indices]) for k, v in stats.max_block_length.items()}),
        mean_block_length=MappingProxyType({k: _readonly(v[indices]) for k, v in stats.mean_block_length.items()}),
    )


def sequence_stats(chains, *, progress=None) -> SequenceStats:
    root = getattr(chains, "_analysis_root", chains)
    indices = getattr(chains, "_root_indices", None)
    cached = getattr(root, "_sequence_stats_cache", None)
    if cached is not None:
        if root is chains or indices is None:
            return cached
        return _slice_sequence_stats(cached, indices)
    # Always compute the expensive statistics on the unfiltered root population.
    # Filtered views then receive a cheap read-only slice of the shared result.
    target = root
    seqs = _sequence_rows(target)
    names = tuple(target.composition.names)
    n = len(seqs)
    trans = np.zeros(n, dtype=np.uint64)
    frac = np.zeros(n, dtype=float)
    block_count = {name: np.zeros(n, dtype=np.uint64) for name in names}
    max_block = {name: np.zeros(n, dtype=np.uint64) for name in names}
    total_block = {name: np.zeros(n, dtype=np.uint64) for name in names}
    total = sum(len(s) for s in seqs)
    enabled = (bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000))
    if enabled:
        print(f"[sequence_stats] started")
    state = {}; done = 0
    for i, seq in enumerate(seqs):
        if not seq:
            continue
        current = seq[0]; run = 1
        for symbol in seq[1:]:
            if symbol == current:
                run += 1
            else:
                trans[i] += 1
                if current in block_count:
                    block_count[current][i] += 1
                    total_block[current][i] += run
                    max_block[current][i] = max(max_block[current][i], run)
                current = symbol; run = 1
        if current in block_count:
            block_count[current][i] += 1
            total_block[current][i] += run
            max_block[current][i] = max(max_block[current][i], run)
        frac[i] = trans[i] / (len(seq) - 1) if len(seq) > 1 else 0.0
        done += len(seq); _progress("sequence_stats", enabled, done, total, state)
    mean_block = {}
    for name in names:
        mean_block[name] = np.divide(total_block[name], block_count[name],
                                     out=np.zeros(n, dtype=float), where=block_count[name] > 0)
    result = SequenceStats(
        names=names,
        transition_count=_readonly(trans), transition_fraction=_readonly(frac),
        block_count=MappingProxyType({k: _readonly(v) for k, v in block_count.items()}),
        max_block_length=MappingProxyType({k: _readonly(v) for k, v in max_block.items()}),
        mean_block_length=MappingProxyType({k: _readonly(v) for k, v in mean_block.items()}),
    )
    try:
        root._sequence_stats_cache = result
    except Exception:
        pass
    if enabled:
        print("[sequence_stats] done")
    if target is not chains and indices is not None:
        return _slice_sequence_stats(result, indices)
    return result


def block_lengths(chains, monomer=None, *, progress=None) -> BlockLengthDistribution:
    seqs = _sequence_rows(chains)
    weights = np.asarray(chains.count, dtype=np.float64)
    counts = {}
    total = sum(len(s) for s in seqs); enabled = (bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000))
    if enabled: print("[block_lengths] started")
    state={}; done=0
    for seq, weight in zip(seqs, weights):
        if seq:
            current=seq[0]; run=1
            for symbol in seq[1:]:
                if symbol == current: run += 1
                else:
                    if monomer is None or current == monomer:
                        counts[run] = counts.get(run, 0.0) + weight
                    current=symbol; run=1
            if monomer is None or current == monomer:
                counts[run] = counts.get(run, 0.0) + weight
        done += len(seq); _progress("block_lengths", enabled, done, total, state)
    lengths=np.asarray(sorted(counts), dtype=np.uint64)
    values=np.asarray([counts[x] for x in lengths], dtype=float)
    fractions=np.divide(values, values.sum(), out=np.zeros_like(values), where=values.sum() > 0)
    if enabled: print("[block_lengths] done")
    return BlockLengthDistribution(monomer, _readonly(lengths), _readonly(values), _readonly(fractions))


def transition_matrix(chains, *, normalize=None, progress=None) -> TransitionMatrix:
    if normalize not in {None, "row", "all"}:
        raise ValueError("normalize must be None, 'row' or 'all'")
    seqs=_sequence_rows(chains); names=tuple(chains.composition.names); index={name:i for i,name in enumerate(names)}
    values=np.zeros((len(names), len(names)), dtype=float); weights=np.asarray(chains.count,dtype=float)
    total=sum(len(s) for s in seqs); enabled=(bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000))
    if enabled: print("[transition_matrix] started")
    state={}; done=0
    for seq, weight in zip(seqs, weights):
        for a,b in zip(seq[:-1],seq[1:]): values[index[a],index[b]] += weight
        done += len(seq); _progress("transition_matrix", enabled, done, total, state)
    if normalize == "row":
        den=values.sum(axis=1,keepdims=True); values=np.divide(values,den,out=np.zeros_like(values),where=den>0)
    elif normalize == "all":
        den=values.sum(); values=values/den if den else values
    if enabled: print("[transition_matrix] done")
    return TransitionMatrix(names, _readonly(values), normalize)


def microstructure_by_dp(chains, statistic, *, monomer=None, bins=None, progress=None) -> MicrostructureByDP:
    stats=sequence_stats(chains, progress=progress)
    if statistic == "transition_count": values=np.asarray(stats.transition_count,dtype=float)
    elif statistic == "transition_fraction": values=np.asarray(stats.transition_fraction,dtype=float)
    elif statistic in {"block_count","max_block_length","mean_block_length"}:
        if monomer is None: raise ValueError(f"monomer is required for {statistic}")
        values=np.asarray(getattr(stats, statistic)[monomer],dtype=float)
    else: raise ValueError("unsupported statistic")
    dp=np.asarray(chains.dp,dtype=float); weights=np.asarray(chains.count,dtype=float); edges=_bin_edges(dp,bins,integer=True)
    n=len(edges)-1; idx=np.searchsorted(edges,dp,side="right")-1; idx[dp==edges[-1]]=n-1
    valid=(idx>=0)&(idx<n)&np.isfinite(values)&(weights>0)
    rc=np.bincount(idx[valid],minlength=n).astype(np.uint64); cc=np.bincount(idx[valid],weights=weights[valid],minlength=n)
    mean=np.full(n,np.nan); med=np.full(n,np.nan); q25=np.full(n,np.nan); q75=np.full(n,np.nan)
    for b in range(n):
        use=valid&(idx==b)
        if np.any(use):
            vv=values[use]; ww=weights[use]; mean[b]=np.average(vv,weights=ww)
            order=np.argsort(vv,kind="stable"); vv=vv[order]; ww=ww[order]; cs=np.cumsum(ww)
            for q,out in ((.25,q25),(.5,med),(.75,q75)):
                out[b]=vv[min(np.searchsorted(cs,q*cs[-1],side="left"),len(vv)-1)]
    return MicrostructureByDP(statistic,monomer,_readonly(edges[:-1]),_readonly(edges[1:]),
        _readonly((edges[:-1]+edges[1:])/2),_readonly(rc),_readonly(cc),_readonly(mean),_readonly(med),_readonly(q25),_readonly(q75))

@dataclass(frozen=True)
class MotifCounts:
    motif: tuple[str, ...]
    occurrence_count: np.ndarray
    normalized_frequency: np.ndarray


@dataclass(frozen=True)
class NGramDistribution:
    n: int
    motifs: tuple[tuple[str, ...], ...]
    count: np.ndarray
    fraction: np.ndarray

    def plot(self, path=None, *, top=None, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        order = np.argsort(self.count)[::-1]
        if top is not None:
            order = order[:int(top)]
        labels = ["|".join(self.motifs[i]) for i in order]
        values = self.count[order]
        ax.bar(np.arange(len(order)), values)
        ax.set_xticks(np.arange(len(order)), labels, rotation=45, ha="right")
        ax.set_xlabel(f"{self.n}-gram")
        ax.set_ylabel("count")
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class PositionProfile:
    names: tuple[str, ...]
    position_left: np.ndarray
    position_right: np.ndarray
    position_center: np.ndarray
    fraction: Mapping[str, np.ndarray]
    mer_count: np.ndarray

    def plot(self, path=None, *, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        for name in self.names:
            ax.plot(self.position_center, self.fraction[name], label=name)
        ax.set_xlabel("relative chain position")
        ax.set_ylabel("monomer fraction")
        ax.set_xlim(0.0, 1.0); ax.set_ylim(0.0, 1.0)
        ax.legend()
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class MicrostructureMap:
    statistic: str
    monomer: str | None
    x_edges: np.ndarray
    y_edges: np.ndarray
    values: np.ndarray

    def plot(self, path=None, *, log=False, style="screen", ax=None, span=None, dpi=300, title=None):
        fig, ax = _axes(ax, style=style, span=span, title=title)
        data = np.log10(self.values + 1.0) if log else self.values
        mesh = ax.pcolormesh(self.x_edges, self.y_edges, data.T, shading="auto")
        ax.set_xlabel("DP")
        label = self.statistic if self.monomer is None else f"{self.statistic} ({self.monomer})"
        ax.set_ylabel(label.replace("_", " "))
        fig.colorbar(mesh, ax=ax, label="log10(chain count + 1)" if log else "chain count")
        fig.tight_layout(); _save(fig, path, dpi)
        return ax


def motif_counts(chains, motif, *, progress=None) -> MotifCounts:
    motif = chains._motif_tokens(motif)
    seqs = _sequence_rows(chains)
    n = len(seqs)
    counts = np.zeros(n, dtype=np.uint64)
    freq = np.zeros(n, dtype=float)
    total = sum(len(s) for s in seqs)
    enabled = bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000)
    if enabled: print("[motif_counts] started")
    state = {}; done = 0; k = len(motif)
    for i, seq in enumerate(seqs):
        possible = max(0, len(seq) - k + 1)
        if possible:
            counts[i] = sum(tuple(seq[j:j+k]) == motif for j in range(possible))
            freq[i] = counts[i] / possible
        done += len(seq); _progress("motif_counts", enabled, done, total, state)
    if enabled: print("[motif_counts] done")
    return MotifCounts(motif, _readonly(counts), _readonly(freq))


def ngrams(chains, n=4, *, min_count=1, progress=None) -> NGramDistribution:
    n = int(n)
    if n < 1:
        raise ValueError("n must be >= 1")
    seqs = _sequence_rows(chains)
    weights = np.asarray(chains.count, dtype=float)
    counts = {}
    total = sum(len(s) for s in seqs)
    enabled = bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000)
    if enabled: print("[ngrams] started")
    state = {}; done = 0
    for seq, weight in zip(seqs, weights):
        for i in range(max(0, len(seq) - n + 1)):
            key = tuple(seq[i:i+n])
            counts[key] = counts.get(key, 0.0) + weight
        done += len(seq); _progress("ngrams", enabled, done, total, state)
    items = sorted((key, value) for key, value in counts.items() if value >= min_count)
    motifs = tuple(key for key, _ in items)
    values = np.asarray([value for _, value in items], dtype=float)
    fractions = np.divide(values, values.sum(), out=np.zeros_like(values), where=values.sum() > 0)
    if enabled: print("[ngrams] done")
    return NGramDistribution(n, motifs, _readonly(values), _readonly(fractions))


def position_profile(chains, *, bins=20, progress=None) -> PositionProfile:
    bins = int(bins)
    if bins < 1:
        raise ValueError("bins must be >= 1")
    seqs = _sequence_rows(chains)
    names = tuple(chains.composition.names)
    index = {name: i for i, name in enumerate(names)}
    values = np.zeros((bins, len(names)), dtype=float)
    totals = np.zeros(bins, dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    total = sum(len(s) for s in seqs)
    enabled = bool(progress) if progress is not None else (bool(options.progress) and total >= 1_000_000)
    if enabled: print("[position_profile] started")
    state = {}; done = 0
    for seq, weight in zip(seqs, weights):
        length = len(seq)
        if length:
            positions = (np.arange(length, dtype=float) + 0.5) / length
            idx = np.minimum((positions * bins).astype(int), bins - 1)
            for b, symbol in zip(idx, seq):
                values[b, index[symbol]] += weight
                totals[b] += weight
        done += length; _progress("position_profile", enabled, done, total, state)
    fractions = np.divide(values, totals[:, None], out=np.zeros_like(values), where=totals[:, None] > 0)
    edges = np.linspace(0.0, 1.0, bins + 1)
    result = PositionProfile(
        names, _readonly(edges[:-1]), _readonly(edges[1:]), _readonly((edges[:-1] + edges[1:]) / 2),
        MappingProxyType({name: _readonly(fractions[:, i]) for i, name in enumerate(names)}),
        _readonly(totals),
    )
    if enabled: print("[position_profile] done")
    return result


def microstructure_map(chains, statistic, *, monomer=None, dp_bins=None, value_bins=None, progress=None) -> MicrostructureMap:
    stats = sequence_stats(chains, progress=progress)
    if statistic == "transition_count":
        values = np.asarray(stats.transition_count, dtype=float)
    elif statistic == "transition_fraction":
        values = np.asarray(stats.transition_fraction, dtype=float)
    elif statistic in {"block_count", "max_block_length", "mean_block_length"}:
        if monomer is None:
            raise ValueError(f"monomer is required for {statistic}")
        values = np.asarray(getattr(stats, statistic)[monomer], dtype=float)
    else:
        raise ValueError("unsupported statistic")
    dp = np.asarray(chains.dp, dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    x_edges = _bin_edges(dp, dp_bins, integer=True)
    y_edges = _bin_edges(values, value_bins, integer=statistic != "transition_fraction")
    valid = np.isfinite(dp) & np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    hist, _, _ = np.histogram2d(dp[valid], values[valid], bins=(x_edges, y_edges), weights=weights[valid])
    return MicrostructureMap(statistic, monomer, _readonly(x_edges), _readonly(y_edges), _readonly(hist))
