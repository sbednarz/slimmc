# copo_stats.nim
# Chain summaries, moments, composition and memory estimates.

import math, tables
import copo_types
import copo_sequence
import ../../common/safe_numeric

proc poolAcceptsChain*(m: Model; poolId: int; c: LiveChain): bool =
  ## True when a live chain is compatible with the terminal/penultimate
  ## metadata declared or inferred for an active pool. A penultimate value
  ## of -1 means "not constrained" and -2 means "mixed penultimate".
  ## Lives here (not in copo_kmc.nim) so both the KMC loop's own channel
  ## eligibility (eligiblePoolIndices in copo_kmc.nim) and the mandatory
  ## validation's pool-consistency check (copo_io.nim, item 3) call the
  ## exact same function -- they can never silently define "belongs to
  ## this pool" two different ways.
  if poolId < 0 or poolId >= m.pools.len: return false
  if m.pools[poolId].kind != pkActive: return false
  if poolId >= m.poolTerminalMer.len: return true
  let terminal = m.poolTerminalMer[poolId]
  if terminal >= 0 and c.last != terminal: return false
  if poolId < m.poolPenultimateMer.len:
    let penultimate = m.poolPenultimateMer[poolId]
    if penultimate >= 0 and c.prev != penultimate: return false
  result = true


proc rawEndgroupMass*(m: Model; name: string): float =
  ## Return a declared end-group contribution independent of the selected mass policy.
  ## Missing end groups are intentionally zero for compatibility.
  for eg in m.endgroups:
    if eg.name == name:
      return eg.mw
  return 0.0

proc endgroupMass*(m: Model; name: string): float =
  ## Return the end-group contribution for the selected model mass policy.
  if m.mass_model != mmWithEndgroups:
    return 0.0
  return m.rawEndgroupMass(name)

proc repeatUnitMass*(m: Model; nMer: seq[int32]): float =
  for i, n in nMer:
    result += float(n) * m.monomers[i].mw

proc recomputeChainMass*(m: Model; nMer: seq[int32]; left_end, right_end: string): float =
  result = m.repeatUnitMass(nMer)
  if m.mass_model == mmWithEndgroups:
    result += m.rawEndgroupMass(left_end) + m.rawEndgroupMass(right_end)

proc countFromConc*(c: float; V: float): int64 =
  result = checkedCountFromConc(c, NA, V, "concentration")

proc concFromCount*(n: int64; V: float): float =
  if V <= 0.0: return 0.0
  float(n) / (NA * V)

proc makeLiveChain*(id: int64; left_end: string; monomerId: int; monomerMw: float;
                    nMonomers: int): LiveChain =
  result.id = id
  result.left_end = left_end
  result.mers = initSequence(monomerId)
  result.right_end = "ACTIVE"
  result.nMer = newSeq[int32](nMonomers)
  result.nMer[monomerId] = 1
  result.dp = 1
  result.mass = monomerMw
  result.last = monomerId
  result.prev = -1
  result.formedBy = fbInit

proc makeLiveChain*(m: Model; id: int64; left_end: string; monomerId: int): LiveChain =
  ## Model-aware constructor. In repeat_units mode this is identical to the
  ## alternate constructor. In with_end_groups mode it adds left_end + ACTIVE.
  result = makeLiveChain(id, left_end, monomerId, m.monomers[monomerId].mw, m.monomers.len)
  if m.mass_model == mmWithEndgroups:
    result.mass += m.endgroupMass(left_end) + m.endgroupMass("ACTIVE")

proc refreshTerminals*(c: var LiveChain) =
  c.last = c.mers.last()
  c.prev = c.mers.prev()

proc pushMonomer*(c: var LiveChain; monomerId: int; monomerMw: float) =
  c.mers.push(monomerId)
  c.dp = checkedAddInt32(c.dp, 1'i32, "chain DP")
  c.nMer[monomerId] = checkedAddInt32(c.nMer[monomerId], 1'i32, "chain composition")
  c.mass += monomerMw
  c.refreshTerminals()

proc popTerminalMonomer*(c: var LiveChain; monomerMw: float): int =
  ## Remove and return the terminal mer. Caller is responsible for returning
  ## the corresponding free monomer count to the simulation state.
  result = c.mers.pop()
  c.dp -= 1
  c.nMer[result] -= 1
  c.mass -= monomerMw
  c.refreshTerminals()

proc makeDeadSummary*(c: LiveChain; right_end: string; formedBy: FormationKind;
                      monomerNames: seq[string]; storeSequence: bool): DeadSummary =
  let nMon = c.nMer.len
  result.left_end = c.left_end
  result.right_end = right_end
  result.dp = c.dp
  result.nMer = c.nMer
  result.mass = c.mass
  result.formedBy = formedBy
  result.firstMer = c.mers.first()
  result.penultimateMer = c.mers.prev()
  result.lastMer = c.mers.last()
  result.dyads = dyadCounts(c.mers, nMon)
  result.triads = triadCounts(c.mers, nMon)
  result.blockCounts = blockCounts(c.mers, nMon)
  result.sequenceStored = storeSequence
  if result.sequenceStored:
    result.sequenceText = c.mers.toText(monomerNames)
  else:
    result.sequenceText = ""
  result.count = 1

proc makeDeadSummary*(m: Model; c: LiveChain; right_end: string; formedBy: FormationKind;
                      monomerNames: seq[string]; storeSequence: bool): DeadSummary =
  ## Convert a live chain into a dead summary and replace the active right-end
  ## mass by the emitted right-end mass when with_end_groups is active.
  result = makeDeadSummary(c, right_end, formedBy, monomerNames, storeSequence)
  if m.mass_model == mmWithEndgroups:
    result.mass = c.mass - m.endgroupMass(c.right_end) + m.endgroupMass(right_end)

proc dyadEndpointBalanceOk*(dyads: openArray[int64]; nMonomers, firstMer,
                            lastMer: int): bool =
  ## For a binary A/B linear sequence, the exact path identity is
  ## N_AB - N_BA = +1/-1/0 according to its endpoints.
  ## This remains checkable on compact DeadSummary records after sequenceText
  ## is dropped, because dyads and endpoints are retained explicitly.
  if nMonomers != 2 or dyads.len != 4:
    return false
  if firstMer < 0 or firstMer >= nMonomers or lastMer < 0 or lastMer >= nMonomers:
    return false
  for a in 0 ..< nMonomers:
    for b in 0 ..< nMonomers:
      if a == b: continue
      let expected =
        if firstMer == a and lastMer == b: 1'i64
        elif firstMer == b and lastMer == a: -1'i64
        else: 0'i64
      if dyads[a * nMonomers + b] - dyads[b * nMonomers + a] != expected:
        return false
  result = true

proc combineLiveToDead*(c1, c2: LiveChain; monomerNames: seq[string];
                        storeSequence: bool): DeadSummary =
  var s = c1.mers
  s.appendReverse(c2.mers)
  var tmp = LiveChain(
    id: -1,
    left_end: c1.left_end,
    mers: s,
    right_end: c2.left_end,
    nMer: newSeq[int32](c1.nMer.len),
    dp: checkedAddInt32(c1.dp, c2.dp, "combined chain DP"),
    mass: c1.mass + c2.mass,
    last: s.last(),
    prev: s.prev(),
    formedBy: fbTermC
  )
  for i in 0 ..< tmp.nMer.len:
    tmp.nMer[i] = checkedAddInt32(c1.nMer[i], c2.nMer[i], "combined chain composition")
  result = makeDeadSummary(tmp, c2.left_end, fbTermC, monomerNames, storeSequence)

proc combineLiveToDead*(m: Model; c1, c2: LiveChain; monomerNames: seq[string];
                        storeSequence: bool): DeadSummary =
  ## Combination joins two live radicals. In with_end_groups mode the two
  ## ACTIVE right-end contributions disappear; the second chain left-end
  ## becomes the right end and its mass contribution is already present.
  result = combineLiveToDead(c1, c2, monomerNames, storeSequence)
  if m.mass_model == mmWithEndgroups:
    result.mass = c1.mass + c2.mass - m.endgroupMass(c1.right_end) - m.endgroupMass(c2.right_end)

proc chainMassForMode*(m: Model; nMer: seq[int32]; left_end, right_end: string; withEndgroups: bool): float =
  result = m.repeatUnitMass(nMer)
  if withEndgroups:
    result += m.rawEndgroupMass(left_end) + m.rawEndgroupMass(right_end)

proc addMoment(ms: var MomentStats; dp: int32; mass: float; count: int64 = 1) =
  if count <= 0: return
  let c = float(count)
  ms.nChains += count
  ms.sumDP += float(dp) * c
  ms.sumDP2 += float(dp) * float(dp) * c
  ms.sumMass += mass * c
  ms.sumMass2 += mass * mass * c
  ms.sumMass3 += mass * mass * mass * c

proc finishMoments(ms: var MomentStats) =
  if ms.nChains > 0:
    ms.dpn = ms.sumDP / float(ms.nChains)
    ms.mn = ms.sumMass / float(ms.nChains)
  if ms.sumDP > 0.0:
    ms.dpw = ms.sumDP2 / ms.sumDP
  if ms.sumMass > 0.0:
    ms.mw = ms.sumMass2 / ms.sumMass
  if ms.sumMass2 > 0.0:
    ms.mz = ms.sumMass3 / ms.sumMass2
  if ms.mn > 0.0:
    ms.pdi = ms.mw / ms.mn

proc structuralMassComplete*(m: Model; s: State): bool
  ## Forward declaration: implemented below, used by computeMomentsForMass
  ## above its own definition to decide whether to NaN-poison endgroup-
  ## inclusive moments.

proc computeMomentsForMass*(m: Model; s: State; withEndgroups: bool; population: string = "all"): MomentStats =
  ## Compute chain moments for a selected population:
  ## - "all": live + dead chains
  ## - "live": active chains only
  ## - "dead": dead summaries only
  ## This makes the output contract explicit before mapper work.
  if population != "dead":
    for p in s.livePools:
      for c in p:
        result.addMoment(c.dp, m.chainMassForMode(c.nMer, c.left_end, c.right_end, withEndgroups))
  if population != "live":
    for d in s.deadChains:
      result.addMoment(d.dp, m.chainMassForMode(d.nMer, d.left_end, d.right_end, withEndgroups), d.count)
  result.finishMoments()
  # Parity fix with classic slimmc: when at least one observed end label has
  # no declared endgroup mass, its true contribution is unknown -- not zero.
  # Classic slimmc poisons its endgroup-inclusive moments with NaN in this
  # situation so the incompleteness is visible directly in the same table a
  # reader is already looking at, rather than only in a side-channel audit
  # file. Only the endgroup-inclusive moments are poisoned; the pure
  # repeat-unit moments (withEndgroups=false) never involve endgroup mass and
  # stay valid regardless of declaration completeness.
  if withEndgroups and not structuralMassComplete(m, s):
    result.mn = NaN
    result.mw = NaN
    result.mz = NaN
    result.pdi = NaN

proc computeMoments*(m: Model; s: State): MomentStats =
  result = computeMomentsForMass(m, s, m.mass_model == mmWithEndgroups, "all")

proc declaredEndgroup*(m: Model; name: string): bool =
  m.endgroupByName.hasKey(name)

proc structuralMassComplete*(m: Model; s: State): bool =
  ## True if every end label currently present on live/dead chains has a
  ## declared endgroup mass. Missing endgroups are still interpreted as zero
  ## in chainMassForMode when the repeat_units mass model is selected.
  ## computeMomentsForMass uses
  ## this flag to NaN-poison its endgroup-inclusive moments instead of
  ## silently reporting a mass total that is known to be wrong.
  for p in s.livePools:
    for c in p:
      if not m.declaredEndgroup(c.left_end): return false
      if not m.declaredEndgroup(c.right_end): return false
  for d in s.deadChains:
    if not m.declaredEndgroup(d.left_end): return false
    if not m.declaredEndgroup(d.right_end): return false
  result = true

proc missingEndgroups*(m: Model; s: State): seq[string] =
  var seen: seq[string] = @[]
  proc addMissing(name: string) =
    if name.len == 0: return
    if m.declaredEndgroup(name): return
    for old in seen:
      if old == name: return
    seen.add name
  for p in s.livePools:
    for c in p:
      addMissing(c.left_end)
      addMissing(c.right_end)
  for d in s.deadChains:
    addMissing(d.left_end)
    addMissing(d.right_end)
  result = seen

proc totalLiveChains*(s: State): int64 =
  for p in s.livePools:
    result += int64(p.len)

proc totalLiveMers*(s: State): int64 =
  for p in s.livePools:
    for c in p:
      result += int64(c.mers.len)

proc totalDeadMersStored*(s: State): int64 =
  for d in s.deadChains:
    if d.sequenceStored:
      result += int64(d.dp) * d.count

proc estimateMemory*(m: Model; s: State): MemoryEstimate =
  result.liveChains = s.totalLiveChains()
  result.deadSummaries = int64(s.deadChains.len)
  result.storedLiveMers = s.totalLiveMers()
  result.storedDeadMers = s.totalDeadMersStored()

  # Conservative, transparent estimates. They are not OS RSS.
  result.liveSeqBytes = result.storedLiveMers       # v0.1: 1 byte / mer
  result.liveObjectBytes = result.liveChains * 96
  var deadBlockCounters = 0'i64
  for d in s.deadChains:
    deadBlockCounters += int64(d.blockCounts.len)
  result.deadSummaryBytes = result.deadSummaries * int64(128 + 16 * m.monomers.len) + deadBlockCounters * 24
  result.totalBytes = result.liveSeqBytes + result.liveObjectBytes + result.deadSummaryBytes

proc addComposition*(acc: var seq[int64]; nMer: seq[int32]; count: int64 = 1) =
  if acc.len < nMer.len: acc.setLen(nMer.len)
  for i in 0 ..< nMer.len:
    acc[i] += int64(nMer[i]) * count

proc polymerCompositionCounts*(m: Model; s: State): seq[int64] =
  result = newSeq[int64](m.monomers.len)
  for p in s.livePools:
    for c in p:
      result.addComposition(c.nMer)
  for d in s.deadChains:
    result.addComposition(d.nMer, d.count)

proc globalDyads*(m: Model; s: State): seq[int64] =
  let n = m.monomers.len
  result = newSeq[int64](n * n)
  for p in s.livePools:
    for c in p:
      let dc = dyadCounts(c.mers, n)
      for i in 0 ..< result.len: result[i] += dc[i]
  for d in s.deadChains:
    for i in 0 ..< result.len: result[i] += d.dyads[i] * d.count

proc globalTriads*(m: Model; s: State): seq[int64] =
  let n = m.monomers.len
  result = newSeq[int64](n * n * n)
  for p in s.livePools:
    for c in p:
      let tc = triadCounts(c.mers, n)
      for i in 0 ..< result.len: result[i] += tc[i]
  for d in s.deadChains:
    for i in 0 ..< result.len: result[i] += d.triads[i] * d.count


proc addBlockCounts(acc: var seq[BlockCount]; counts: seq[BlockCount]; multiplier: int64 = 1) =
  if multiplier <= 0: return
  for bc in counts:
    var found = false
    for i in 0 ..< acc.len:
      if acc[i].monomerId == bc.monomerId and acc[i].length == bc.length:
        acc[i].count += bc.count * multiplier
        found = true
        break
    if not found:
      acc.add BlockCount(monomerId: bc.monomerId, length: bc.length, count: bc.count * multiplier)

proc globalBlockCounts*(m: Model; s: State): seq[BlockCount] =
  ## Aggregate exact contiguous-run histograms from live full sequences and
  ## dead compact summaries.  Dead summaries preserve blockCounts even when
  ## sequenceText was dropped.
  let n = m.monomers.len
  for p in s.livePools:
    for c in p:
      addBlockCounts(result, blockCounts(c.mers, n))
  for d in s.deadChains:
    addBlockCounts(result, d.blockCounts, d.count)
