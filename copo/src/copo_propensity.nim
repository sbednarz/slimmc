# Exact copo propensity and eligibility helpers shared by the SSA and Storage snapshots.
import math
import copo_types
import copo_stats

proc eligiblePoolIndices*(m: Model; s: State; poolId: int): seq[int] =
  ## Return indices of chains that really belong to the active-pool class.
  ## This hardens penultimate channels against stale hidden PA/PB assumptions.
  if poolId < 0 or poolId >= s.livePools.len: return
  for i, c in s.livePools[poolId]:
    if m.poolAcceptsChain(poolId, c):
      result.add i

proc eligiblePropIndices*(m: Model; s: State; poolId: int): seq[int] =
  if poolId < 0 or poolId >= s.livePools.len: return
  for i, c in s.livePools[poolId]:
    if m.poolAcceptsChain(poolId, c) and int64(c.dp) < m.dp_max:
      result.add i

proc countersCurrent(s: State; poolId: int): bool =
  poolId >= 0 and
  poolId < s.livePools.len and
  poolId < s.poolEligibleCounts.len and
  poolId < s.poolPropagatableCounts.len and
  poolId < s.poolEligibilityTrackedLen.len and
  s.poolEligibilityTrackedLen[poolId] == int64(s.livePools[poolId].len)

proc eligiblePropCount*(m: Model; s: State; poolId: int): int =
  if poolId < 0 or poolId >= s.livePools.len: return 0
  if s.countersCurrent(poolId):
    return int(s.poolPropagatableCounts[poolId])
  for c in s.livePools[poolId]:
    if m.poolAcceptsChain(poolId, c) and int64(c.dp) < m.dp_max:
      inc result

proc eligiblePoolCount*(m: Model; s: State; poolId: int): int =
  if poolId < 0 or poolId >= s.livePools.len: return 0
  if s.countersCurrent(poolId):
    return int(s.poolEligibleCounts[poolId])
  for c in s.livePools[poolId]:
    if m.poolAcceptsChain(poolId, c):
      inc result

proc eligibleDepropCount*(m: Model; s: State; ch: KmcChannel): int =
  if ch.pool1 < 0 or ch.pool1 >= s.livePools.len: return 0
  if ch.poolOut < 0 or ch.poolOut >= m.poolTerminalMer.len: return 0
  let newTerminal = m.poolTerminalMer[ch.poolOut]
  if newTerminal < 0: return 0
  for c in s.livePools[ch.pool1]:
    if c.dp >= 2 and c.last == ch.monomerId and c.prev == newTerminal:
      inc result

proc eligibleDepropIndices*(m: Model; s: State; ch: KmcChannel): seq[int] =
  ## A deprop channel is specific to both the removed terminal mer and the
  ## active pool that receives the shortened chain.
  ## Example: macro deprop PA -> PB + A kBA handles ...B-A* -> ...B* + A.
  if ch.pool1 < 0 or ch.pool1 >= s.livePools.len: return
  if ch.poolOut < 0 or ch.poolOut >= m.poolTerminalMer.len: return
  let newTerminal = m.poolTerminalMer[ch.poolOut]
  if newTerminal < 0: return
  for i, c in s.livePools[ch.pool1]:
    if c.dp >= 2 and c.last == ch.monomerId and c.prev == newTerminal:
      result.add i


proc smallCount(s: State; r: SmallRef): int64 =
  case r.kind
  of skSpecies: result = s.speciesN[r.id]
  of skMonomer: result = s.monomerN[r.id]

proc falling(n: int64; k: int): float =
  if k <= 0: return 1.0
  if n < int64(k): return 0.0
  result = 1.0
  for i in 0 ..< k:
    result *= float(n - int64(i))

proc elementaryPropensity(m: Model; s: State; ch: KmcChannel): float =
  # Items 79/80/81: propensity is driven by k alone, never multiplied by
  # f. f only gates whether products form once a firing is already
  # selected by the SSA -- it does not change how often the reaction is
  # attempted. (Previously this multiplied by ch.efficiency, silently
  # halving the attempt rate for f<1 reactions instead of the correct
  # "always attempt at k, sometimes fail to produce" semantics.)
  let k = m.rateValue(ch.kId)
  if ch.kind == chRxnUni:
    let n = smallCount(s, ch.smallReactants[0])
    result = k * float(n)
  elif ch.kind == chRxnBiSame:
    let r = ch.smallReactants[0]
    let n = smallCount(s, r)
    result = k / (NA * m.V) * falling(n, 2)
  elif ch.kind == chRxnBiDiff:
    let n1 = smallCount(s, ch.smallReactants[0])
    let n2 = smallCount(s, ch.smallReactants[1])
    result = k / (NA * m.V) * float(n1) * float(n2)

proc computePropensities*(m: Model; s: State): seq[float] =
  result = newSeq[float](m.channels.len)
  for i, ch in m.channels:
    let k = m.rateValue(ch.kId)
    case ch.kind
    of chRxnUni, chRxnBiDiff, chRxnBiSame:
      result[i] = elementaryPropensity(m, s, ch)
    of chMacroInit:
      let nr = s.speciesN[ch.speciesId]
      let nm = s.monomerN[ch.monomerId]
      if nr > 0 and nm > 0:
        result[i] = k / (NA * m.V) * float(nr) * float(nm)
    of chMacroProp:
      let np = m.eligiblePropCount(s, ch.pool1)
      let nm = s.monomerN[ch.monomerId]
      if np > 0 and nm > 0:
        result[i] = k / (NA * m.V) * float(np) * float(nm)
    of chMacroTermC, chMacroTermD:
      let n1 = m.eligiblePoolCount(s, ch.pool1)
      let n2 = m.eligiblePoolCount(s, ch.pool2)
      if ch.pool1 == ch.pool2:
        if n1 >= 2:
          result[i] = k / (NA * m.V) * float(n1) * float(n1 - 1)
      else:
        if n1 > 0 and n2 > 0:
          result[i] = k / (NA * m.V) * float(n1) * float(n2)
    of chMacroTermX:
      let np = m.eligiblePoolCount(s, ch.pool1)
      let nx = s.speciesN[ch.speciesId]
      if np > 0 and nx > 0:
        result[i] = k / (NA * m.V) * float(np) * float(nx)
    of chMacroTransfer:
      let np = m.eligiblePoolCount(s, ch.pool1)
      let nx = s.speciesN[ch.speciesId]
      if np > 0 and nx > 0:
        result[i] = k / (NA * m.V) * float(np) * float(nx)
    of chMacroTransferM:
      let np = m.eligiblePoolCount(s, ch.pool1)
      let nm = s.monomerN[ch.monomerId]
      if np > 0 and nm > 0:
        result[i] = k / (NA * m.V) * float(np) * float(nm)
    of chMacroDeprop:
      let nElig =
        if s.countersCurrent(ch.pool1) and s.channelDepropEligibleCounts.len == m.channels.len:
          int(s.channelDepropEligibleCounts[i])
        else:
          m.eligibleDepropCount(s, ch)
      if nElig > 0:
        result[i] = k * float(nElig)

