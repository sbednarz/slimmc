# slimmc_kmc.nim
# SSA/kMC engine for slimmc.

import math, random, times, tables, strutils, os
import slimmc_types
import slimmc_parser
import slimmc_io
import slimmc_storage
import ../../common/results_types
import ../../common/safe_numeric

proc initState*(m: Model): State =
  result.t = 0.0
  result.kmcEvent = 0
  result.actionNo = 0
  result.scheduledActionNo = 0
  result.conditionalActionNo = 0
  result.stateRevision = 0
  result.parameterStateId = 1
  result.snapshotId = 0
  result.snapshots = @[]
  result.anySnapshotWritten = false
  result.chainsWritten = false
  result.channelTraceRowsWritten = 0
  result.channelTraceTruncated = false
  result.n = newSeq[int64](m.species.len)
  result.speciesExternalN = newSeq[int64](m.species.len)
  result.speciesDosedN = newSeq[int64](m.species.len)
  result.pools = newSeq[seq[Chain]](m.pools.len)
  result.poolPropagatableCounts = newSeq[int64](m.pools.len)
  result.poolDepropableCounts = newSeq[int64](m.pools.len)
  result.poolEligibilityTrackedLen = newSeq[int64](m.pools.len)
  result.channelFires = newSeq[int64](m.channels.len)
  result.channelSuccesses = newSeq[int64](m.channels.len)
  result.channelFailures = newSeq[int64](m.channels.len)
  result.observedEg = newSeq[bool](m.egNames.len)
  result.feedRemainders = newSeq[seq[float]](m.feeds.len)
  for fid in 0 ..< m.feeds.len:
    result.feedRemainders[fid] = newSeq[float](m.species.len)

  for i, sp in m.species:
    result.n[i] = countFromConc(sp.c0, m.V)
    result.speciesExternalN[i] = result.n[i]

  if m.monomerId >= 0:
    result.mExpected = result.n[m.monomerId]
    result.mBalance = result.n[m.monomerId]

proc eligibilityCountersCurrent(s: State; poolId: int): bool =
  poolId >= 0 and
  poolId < s.pools.len and
  poolId < s.poolPropagatableCounts.len and
  poolId < s.poolDepropableCounts.len and
  poolId < s.poolEligibilityTrackedLen.len and
  s.poolEligibilityTrackedLen[poolId] == int64(s.pools[poolId].len)

proc countPropagatable(s: State; poolId: int; dpMax: int64): int64 =
  if poolId < 0 or poolId >= s.pools.len: return 0
  if s.eligibilityCountersCurrent(poolId):
    return s.poolPropagatableCounts[poolId]
  for ch in s.pools[poolId]:
    if int64(ch.dp) < dpMax:
      inc result

proc choosePropagatable(s: State; poolId: int; dpMax: int64; rng: var Rand): int =
  let pool = s.pools[poolId]
  let nOk = s.countPropagatable(poolId, dpMax)
  doAssert nOk > 0
  if nOk == int64(pool.len):
    return rng.rand(pool.len - 1)
  let pick = rng.rand(int(nOk) - 1)
  var seen = 0
  for i, ch in pool:
    if int64(ch.dp) < dpMax:
      if seen == pick: return i
      inc seen
  raise newException(AssertionDefect, "propagatable chain not found")

proc countDepropable(s: State; poolId: int): int64 =
  if poolId < 0 or poolId >= s.pools.len: return 0
  if s.eligibilityCountersCurrent(poolId):
    return s.poolDepropableCounts[poolId]
  for ch in s.pools[poolId]:
    if ch.dp > 1:
      inc result

proc chooseDepropable(s: State; poolId: int; rng: var Rand): int =
  let pool = s.pools[poolId]
  let nOk = s.countDepropable(poolId)
  doAssert nOk > 0
  if nOk == int64(pool.len):
    return rng.rand(pool.len - 1)
  let pick = rng.rand(int(nOk) - 1)
  var seen = 0
  for i, ch in pool:
    if ch.dp > 1:
      if seen == pick: return i
      inc seen
  raise newException(AssertionDefect, "depropagatable chain not found")

proc rebuildEligibilityCounters(s: var State; m: Model; poolId: int) =
  if poolId < 0 or poolId >= s.pools.len: return
  if s.poolPropagatableCounts.len != s.pools.len:
    s.poolPropagatableCounts = newSeq[int64](s.pools.len)
  if s.poolDepropableCounts.len != s.pools.len:
    s.poolDepropableCounts = newSeq[int64](s.pools.len)
  if s.poolEligibilityTrackedLen.len != s.pools.len:
    s.poolEligibilityTrackedLen = newSeq[int64](s.pools.len)
  var propCount = 0'i64
  var depropCount = 0'i64
  for ch in s.pools[poolId]:
    if int64(ch.dp) < m.dpMax: inc propCount
    if ch.dp > 1: inc depropCount
  s.poolPropagatableCounts[poolId] = propCount
  s.poolDepropableCounts[poolId] = depropCount
  s.poolEligibilityTrackedLen[poolId] = int64(s.pools[poolId].len)

proc ensureEligibilityCounters(s: var State; m: Model; poolId: int) =
  if not s.eligibilityCountersCurrent(poolId):
    s.rebuildEligibilityCounters(m, poolId)

proc computePropensities*(m: Model; s: State): seq[float] =
  result = newSeq[float](m.channels.len)

  for i, ch in m.channels:
    let k = rateValue(m, ch.kId)

    case ch.kind
    of chRxnUni:
      let a = ch.lhs[0]
      if s.n[a.sp] >= int64(a.stoich):
        result[i] = k * float(s.n[a.sp])

    of chRxnBiDiff:
      let a = ch.lhs[0]
      let b = ch.lhs[1]
      if s.n[a.sp] > 0 and s.n[b.sp] > 0:
        result[i] = k / (NA * m.V) * float(s.n[a.sp]) * float(s.n[b.sp])

    of chRxnBiSame:
      let a = ch.lhs[0]
      let n = s.n[a.sp]
      if n >= 2:
        result[i] = k / (NA * m.V) * float(n) * float(n - 1)

    of chMacroInit:
      let nr = s.n[ch.sp1]
      let nm = s.n[ch.sp2]
      if nr > 0 and nm > 0:
        result[i] = k / (NA * m.V) * float(nr) * float(nm)

    of chMacroProp:
      let np = s.countPropagatable(ch.pool1, m.dpMax)
      let nm = s.n[ch.sp1]
      if np > 0 and nm > 0:
        result[i] = k / (NA * m.V) * float(np) * float(nm)

    of chMacroDeprop:
      let nd = s.countDepropable(ch.pool1)
      if nd > 0:
        result[i] = k * float(nd)

    of chMacroTermC, chMacroTermD:
      let n1 = s.pools[ch.pool1].len
      let n2 = s.pools[ch.pool2].len
      if ch.pool1 == ch.pool2:
        if n1 >= 2:
          result[i] = k / (NA * m.V) * float(n1) * float(n1 - 1)
      elif n1 > 0 and n2 > 0:
        result[i] = k / (NA * m.V) * float(n1) * float(n2)

    of chMacroTermX, chMacroTransferH, chMacroTransferM:
      let np = s.pools[ch.pool1].len
      let ns = s.n[ch.sp1]
      if np > 0 and ns > 0:
        result[i] = k / (NA * m.V) * float(np) * float(ns)

proc takeAt(pool: var seq[Chain]; idx: int): Chain =
  result = pool[idx]
  pool[idx] = pool[^1]
  pool.setLen(pool.len - 1)

proc trackAdd(s: var State; m: Model; poolId: int; c: Chain) =
  s.ensureEligibilityCounters(m, poolId)
  s.pools[poolId].add c
  if int64(c.dp) < m.dpMax: inc s.poolPropagatableCounts[poolId]
  if c.dp > 1: inc s.poolDepropableCounts[poolId]
  inc s.poolEligibilityTrackedLen[poolId]

proc trackRemove(s: var State; m: Model; poolId: int; idx: int): Chain =
  s.ensureEligibilityCounters(m, poolId)
  let c = s.pools[poolId][idx]
  if int64(c.dp) < m.dpMax: dec s.poolPropagatableCounts[poolId]
  if c.dp > 1: dec s.poolDepropableCounts[poolId]
  result = s.pools[poolId].takeAt(idx)
  dec s.poolEligibilityTrackedLen[poolId]

proc trackDpChange(s: var State; m: Model; poolId: int; idx: int; delta: int) =
  s.ensureEligibilityCounters(m, poolId)
  let oldDp = s.pools[poolId][idx].dp
  if int64(oldDp) < m.dpMax: dec s.poolPropagatableCounts[poolId]
  if oldDp > 1: dec s.poolDepropableCounts[poolId]
  s.pools[poolId][idx].dp += delta
  let newDp = s.pools[poolId][idx].dp
  if int64(newDp) < m.dpMax: inc s.poolPropagatableCounts[poolId]
  if newDp > 1: inc s.poolDepropableCounts[poolId]

proc untrackSelected(s: var State; m: Model; poolId: int; idx: int) =
  s.ensureEligibilityCounters(m, poolId)
  let c = s.pools[poolId][idx]
  if int64(c.dp) < m.dpMax: dec s.poolPropagatableCounts[poolId]
  if c.dp > 1: dec s.poolDepropableCounts[poolId]

proc chooseTwoDistinct(n: int; rng: var Rand): tuple[i: int, j: int] =
  let i = rng.rand(n - 1)
  var j = rng.rand(n - 2)
  if j >= i:
    inc j
  result = (i: i, j: j)

proc takeTwo(pool: var seq[Chain]; i0, j0: int): tuple[a: Chain, b: Chain] =
  var i = i0
  var j = j0
  if i == j:
    raise newException(ValueError, "takeTwo got identical indices")
  if i > j:
    result.a = pool[i]
    result.b = pool[j]
    discard pool.takeAt(i)
    discard pool.takeAt(j)
  else:
    result.a = pool[i]
    result.b = pool[j]
    discard pool.takeAt(j)
    discard pool.takeAt(i)

proc speciesEndGroup(m: Model; sp: int): int =
  result = m.egByName[endGroupName(m.species[sp].name)]

proc applyChannel*(m: var Model; s: var State; chId: int; rng: var Rand) =
  let ch = m.channels[chId]

  case ch.kind
  of chRxnUni, chRxnBiDiff, chRxnBiSame:
    for tr in ch.lhs:
      s.n[tr.sp] -= int64(tr.stoich)
      if tr.sp == m.monomerId:
        s.mBalance -= int64(tr.stoich)
    if rng.rand(1.0) <= ch.eff:
      for tr in ch.rhs:
        s.n[tr.sp] += int64(tr.stoich)
        if tr.sp == m.monomerId:
          s.mBalance += int64(tr.stoich)
      s.channelSuccesses[chId] += 1
    else:
      s.channelFailures[chId] += 1

  of chMacroInit:
    s.n[ch.sp1] -= 1
    s.n[ch.sp2] -= 1
    s.trackAdd(m, ch.poolOut, Chain(
      eg1: speciesEndGroup(m, ch.sp1), dp: 1, eg2: m.egActive,
      formedBy: fbInit
    ))

  of chMacroProp:
    s.n[ch.sp1] -= 1
    let idx = choosePropagatable(s, ch.pool1, m.dpMax, rng)
    s.trackDpChange(m, ch.pool1, idx, 1)

  of chMacroDeprop:
    let idx = chooseDepropable(s, ch.pool1, rng)
    s.trackDpChange(m, ch.pool1, idx, -1)
    s.n[ch.sp1] += 1

  of chMacroTermC:
    var c1, c2: Chain
    if ch.pool1 == ch.pool2:
      let ij = chooseTwoDistinct(s.pools[ch.pool1].len, rng)
      s.untrackSelected(m, ch.pool1, ij.i)
      s.untrackSelected(m, ch.pool1, ij.j)
      let pair = takeTwo(s.pools[ch.pool1], ij.i, ij.j)
      s.poolEligibilityTrackedLen[ch.pool1] -= 2
      c1 = pair.a
      c2 = pair.b
    else:
      c1 = s.trackRemove(m, ch.pool1, rng.rand(s.pools[ch.pool1].len - 1))
      c2 = s.trackRemove(m, ch.pool2, rng.rand(s.pools[ch.pool2].len - 1))
    s.trackAdd(m, ch.poolOut, Chain(
      eg1: c1.eg1, dp: c1.dp + c2.dp, eg2: c2.eg1,
      formedBy: fbTermC
    ))

  of chMacroTermD:
    var c1, c2: Chain
    if ch.pool1 == ch.pool2:
      let ij = chooseTwoDistinct(s.pools[ch.pool1].len, rng)
      s.untrackSelected(m, ch.pool1, ij.i)
      s.untrackSelected(m, ch.pool1, ij.j)
      let pair = takeTwo(s.pools[ch.pool1], ij.i, ij.j)
      s.poolEligibilityTrackedLen[ch.pool1] -= 2
      c1 = pair.a
      c2 = pair.b
    else:
      c1 = s.trackRemove(m, ch.pool1, rng.rand(s.pools[ch.pool1].len - 1))
      c2 = s.trackRemove(m, ch.pool2, rng.rand(s.pools[ch.pool2].len - 1))
    if rng.rand(1) == 0:
      s.trackAdd(m, ch.poolOut, Chain(
        eg1: c1.eg1, dp: c1.dp, eg2: m.egH, formedBy: fbTermD_H
      ))
      s.trackAdd(m, ch.poolOut, Chain(
        eg1: c2.eg1, dp: c2.dp, eg2: m.egU, formedBy: fbTermD_U
      ))
    else:
      s.trackAdd(m, ch.poolOut, Chain(
        eg1: c1.eg1, dp: c1.dp, eg2: m.egU, formedBy: fbTermD_U
      ))
      s.trackAdd(m, ch.poolOut, Chain(
        eg1: c2.eg1, dp: c2.dp, eg2: m.egH, formedBy: fbTermD_H
      ))

  of chMacroTermX:
    s.n[ch.sp1] -= 1
    let c = s.trackRemove(m, ch.pool1, rng.rand(s.pools[ch.pool1].len - 1))
    s.trackAdd(m, ch.poolOut, Chain(
      eg1: c.eg1, dp: c.dp, eg2: speciesEndGroup(m, ch.sp1),
      formedBy: fbTermX
    ))

  of chMacroTransferH:
    s.n[ch.sp1] -= 1
    s.n[ch.sp2] += 1
    let c = s.trackRemove(m, ch.pool1, rng.rand(s.pools[ch.pool1].len - 1))
    s.trackAdd(m, ch.poolOut, Chain(
      eg1: c.eg1, dp: c.dp, eg2: m.egH, formedBy: fbTransferH
    ))

  of chMacroTransferM:
    # P_n* + M -> D_n-H + P_1*.
    # The consumed monomer is the first repeat unit of the newborn active chain;
    # this is deliberately one SSA firing, not transfer_h followed by macro init.
    # The eg1 label is <monomer>_tr. DP=1 already includes the first repeat
    # unit; its built-in -1.008 g/mol contribution records the transfer-to-
    # monomer structural correction used by the optional end-group mass model.
    s.n[ch.sp1] -= 1
    let c = s.trackRemove(m, ch.pool1, rng.rand(s.pools[ch.pool1].len - 1))
    s.trackAdd(m, ch.poolOut, Chain(
      eg1: c.eg1, dp: c.dp, eg2: m.egH, formedBy: fbTransferM
    ))
    s.trackAdd(m, ch.pool2, Chain(
      eg1: m.egByName[transferMonomerEndGroupName(m.species[ch.sp1].name)],
      dp: 1,
      eg2: m.egActive,
      formedBy: fbTransferM
    ))

proc nextScheduledActionTime(actions: seq[ScheduledAction]): float =
  result = InfTime
  for e in actions:
    if e.active and e.nextTime < result:
      result = e.nextTime

proc executeAction(
  m: var Model;
  s: var State;
  action: ActionKind;
  args: seq[string];
  lineNo: int;
  wallStart: float;
  opts: RunOptions
): ActionResult =
  result.before = 0.0
  result.after = 0.0
  result.outputWritten = ""

  case action
  of eaPrint:
    result.message = args[0]
    printMarker(m, opts, result.message)
  of eaPrintInfo:
    printProgress(m, s, wallStart, opts)

  of eaSave:
    if args.len != 0:
      fail(lineNo, "save takes no arguments")
    saveSnapshot(m, s, withChains = false, source = "save")
    result.outputWritten = "snapshot"

  of eaSaveChains:
    if args.len != 0:
      fail(lineNo, "save_chains takes no arguments")
    saveSnapshot(m, s, withChains = true, source = "save_chains")
    result.outputWritten = "chains"

  of eaStop:
    if args.len != 0:
      fail(lineNo, "stop takes no arguments")
    result.message = "stop requested"

  of eaSetK:
    if args.len != 2:
      fail(lineNo, "set_k syntax: set_k rate value")
    let kId = m.getRateId(args[0], lineNo)
    result.target = args[0]
    result.requested = args[1]
    result.before = rateValue(m, kId)
    result.after = parseF(args[1], lineNo, "action numeric argument")
    if result.after < 0.0:
      fail(lineNo, "set_k requires rate value >= 0")
    m.rates[kId].kind = rkFixed
    m.rates[kId].kConst = result.after
    result.hasNumeric = true
    result.stateChanged = true

  of eaAddK:
    if args.len != 2:
      fail(lineNo, "add_k syntax: add_k rate increment")
    let kId = m.getRateId(args[0], lineNo)
    result.target = args[0]
    result.requested = args[1]
    result.before = rateValue(m, kId)
    result.after = result.before + parseF(args[1], lineNo, "action numeric argument")
    if result.after < 0.0:
      fail(lineNo, "add_k would make negative rate constant")
    m.rates[kId].kind = rkFixed
    m.rates[kId].kConst = result.after
    result.hasNumeric = true
    result.stateChanged = true

  of eaSetTemp:
    if args.len != 1:
      fail(lineNo, "set_temp syntax: set_temp value")
    result.target = "T"
    result.requested = args[0]
    result.before = m.T
    result.after = parseF(args[0], lineNo, "action numeric argument")
    if result.after <= 0.0:
      fail(lineNo, "set_temp requires temperature > 0")
    m.T = result.after
    result.hasNumeric = true
    result.stateChanged = true

  of eaAddTemp:
    if args.len != 1:
      fail(lineNo, "add_temp syntax: add_temp increment")
    result.target = "T"
    result.requested = args[0]
    result.before = m.T
    result.after = result.before + parseF(args[0], lineNo, "action numeric argument")
    if result.after <= 0.0:
      fail(lineNo, "add_temp would make non-positive temperature")
    m.T = result.after
    result.hasNumeric = true
    result.stateChanged = true

  of eaSetC:
    if args.len != 2:
      fail(lineNo, "set_c syntax: set_c species value")
    let sp = m.getSpeciesId(args[0], lineNo)
    let warning = "WARNING: line " & $lineNo & ": set_c forces the concentration of '" & args[0] & "'. Its physical material balance is invalid after this action."
    stdout.writeLine(warning)
    stdout.flushFile()
    logLine(m, opts, warning)
    result.target = args[0]
    result.requested = args[1]
    result.before = conc(s.n[sp], m.V)
    let requested = parseF(args[1], lineNo, "action numeric argument")
    if requested < 0.0:
      fail(lineNo, "set_c requires concentration >= 0")
    let oldN = s.n[sp]
    let newN = checkedCountFromConc(requested, NA, m.V, "set_c")
    let delta = checkedSubInt64(newN, oldN, "set_c count delta")
    s.n[sp] = newN
    if sp == m.monomerId:
      s.mExpected = checkedAddInt64(s.mExpected, delta, "set_c expected monomer balance")
      s.mBalance = checkedAddInt64(s.mBalance, delta, "set_c monomer balance")
    result.after = conc(s.n[sp], m.V)
    result.hasNumeric = true
    result.stateChanged = true

  of eaAddC:
    if args.len != 2:
      fail(lineNo, "add_c syntax: add_c species increment")
    let sp = m.getSpeciesId(args[0], lineNo)
    result.target = args[0]
    result.requested = args[1]
    result.before = conc(s.n[sp], m.V)
    let increment = parseF(args[1], lineNo, "action numeric argument")
    let dn = checkedCountFromConc(increment, NA, m.V, "add_c", allowNegative = true)
    let newN = checkedAddInt64(s.n[sp], dn, "add_c molecule count")
    if newN < 0:
      fail(lineNo, "add_c would make negative molecule count")
    s.n[sp] = newN
    s.speciesExternalN[sp] = checkedAddInt64(s.speciesExternalN[sp], dn, "add_c external balance")
    if sp == m.monomerId:
      s.mExpected = checkedAddInt64(s.mExpected, dn, "add_c expected monomer balance")
      s.mBalance = checkedAddInt64(s.mBalance, dn, "add_c monomer balance")
    result.after = conc(s.n[sp], m.V)
    result.hasNumeric = true
    result.stateChanged = true


  of eaFeed:
    if args.len notin [2, 3]:
      fail(lineNo, "feed action syntax: feed NAME VOLUME [L|l|mL|ml|ML]; without a unit VOLUME is in L")
    if not m.feedByName.hasKey(args[0]):
      fail(lineNo, "unknown feed: " & args[0])
    let doseValue = parseF(args[1], lineNo, "feed volume")
    var doseMl: float
    if args.len == 2:
      doseMl = doseValue * 1000.0
    else:
      let unit = args[2].toLowerAscii()
      case unit
      of "l": doseMl = doseValue * 1000.0
      of "ml": doseMl = doseValue
      else: fail(lineNo, "feed volume unit must be L or mL (case-insensitive)")
    if doseMl <= 0.0:
      fail(lineNo, "feed volume must be > 0")
    if m.currentVolumeMl <= 0.0:
      fail(lineNo, "feed requires param init_volume > 0")
    let fid = m.feedByName[args[0]]
    let oldV = m.V
    let newPhysical = m.currentVolumeMl + doseMl
    let newV = oldV * newPhysical / m.currentVolumeMl
    let deltaV = newV - oldV
    result.target = args[0]
    result.requested = $doseMl
    result.before = m.currentVolumeMl
    for sp in 0 ..< m.species.len:
      let exact = m.feeds[fid].concentrations[sp] * NA * deltaV + s.feedRemainders[fid][sp]
      if exact < 0.0 or exact > float(high(int64)):
        fail(lineNo, "feed molecule count overflow for " & m.species[sp].name)
      let dn = int64(floor(exact))
      s.feedRemainders[fid][sp] = exact - float(dn)
      s.n[sp] = checkedAddInt64(s.n[sp], dn, "feed molecule count")
      s.speciesExternalN[sp] = checkedAddInt64(s.speciesExternalN[sp], dn, "feed external balance")
      s.speciesDosedN[sp] = checkedAddInt64(s.speciesDosedN[sp], dn, "feed dosed balance")
      if sp == m.monomerId:
        s.mExpected = checkedAddInt64(s.mExpected, dn, "feed expected monomer balance")
        s.mBalance = checkedAddInt64(s.mBalance, dn, "feed monomer balance")
    m.currentVolumeMl = newPhysical
    m.V = newV
    result.after = m.currentVolumeMl
    result.hasNumeric = true
    result.stateChanged = true

  of eaPrintMemory:
    if args.len != 0:
      fail(lineNo, "print_memory takes no arguments")
    printMemory(m, s)
    result.outputWritten = "memory"

proc changesParameterState(action: ActionKind): bool =
  action in {eaSetK, eaAddK, eaSetTemp, eaAddTemp}

proc finalizeStateChange(m: Model; s: var State; action: ActionKind) =
  if not changesParameterState(action):
    return
  inc s.parameterStateId
  captureKineticParameterSetV1(m, s, true, uint64(max(0'i64, s.actionNo - 1)))
  writeParameterState(m, s, s.actionNo)
  saveSnapshot(m, s, withChains = false, source = "parameter_change")

proc conditionalValue(m: Model; s: State; a: AtomicCondition): float =
  case a.observable
  of woConversion:
    result = monomerConversion(m, s)
  of woSpeciesConc:
    result = conc(s.n[a.speciesId], m.V)

proc conditionalIsTrue(value: float; a: AtomicCondition): bool =
  case a.comparison
  of coGreater: result = value > a.threshold
  of coLess: result = value < a.threshold

proc conditionText(m: Model; e: ConditionalAction): string =
  var parts: seq[string] = @[]
  for a in e.conditions:
    parts.add conditionalObservableName(m, a) & " " & comparisonName(a.comparison) & " " & num(a.threshold)
  result = parts.join(" and ")

proc conditionValues(values: seq[float]): string =
  var parts: seq[string] = @[]
  for v in values: parts.add num(v)
  result = parts.join(" and ")

proc checkConditionalActions(
  m: var Model;
  s: var State;
  wallStart: float;
  opts: RunOptions;
  checkSource: string
) =
  ## Each line fires once. All atoms on one line are evaluated on the same
  ## state. A stop request is honored only after this complete cascade scan.
  var fired = true
  while fired:
    fired = false
    for i in 0 ..< m.conditionalActions.len:
      if not m.conditionalActions[i].active:
        continue
      let e = m.conditionalActions[i]
      var values: seq[float] = @[]
      var matched = true
      for a in e.conditions:
        let value = conditionalValue(m, s, a)
        values.add value
        if not conditionalIsTrue(value, a): matched = false
      if not matched:
        continue

      m.conditionalActions[i].active = false
      inc s.actionNo
      inc s.conditionalActionNo
      let actionResult = executeAction(m, s, e.action, e.args, e.lineNo, wallStart, opts)
      if actionResult.stateChanged:
        inc s.stateRevision
        finalizeStateChange(m, s, e.action)
      if e.action == eaStop:
        s.stopRequested = true
        s.stopLineNo = e.lineNo
        s.stopCheckSource = checkSource
        s.stopConditions = e.conditions
        s.stopActualValues = values

      writeActionTraceRow(
        m, s, e.lineNo, "when", "", checkSource,
        conditionText(m, e), (if e.conditions.len > 1: "and" else: comparisonName(e.conditions[0].comparison)),
        (if e.conditions.len > 1: conditionText(m, e) else: num(e.conditions[0].threshold)),
        conditionValues(values), e.action, actionResult, e.conditions, values
      )

      if opts.debug:
        debugCheckState(m, s, "when line=" & $e.lineNo & " source=" & checkSource)
      fired = true

proc executeScheduledAction(
  m: var Model;
  s: var State;
  e: ScheduledAction;
  wallStart: float;
  opts: RunOptions
) =
  inc s.actionNo
  inc s.scheduledActionNo
  let actionResult = executeAction(m, s, e.action, e.args, e.lineNo, wallStart, opts)
  if actionResult.stateChanged:
    inc s.stateRevision
    finalizeStateChange(m, s, e.action)

  writeActionTraceRow(
    m, s, e.lineNo, (if e.repeat: "every" else: "at"), num(e.nextTime),
    "time", "", "", "", "", e.action, actionResult
  )

  if opts.debug:
    debugCheckState(m, s, actionKindName(e.action) & " line=" & $e.lineNo)

  if actionResult.stateChanged and m.conditionalActions.len > 0:
    checkConditionalActions(m, s, wallStart, opts, "scheduled_action")

proc processDueActions(m: var Model; s: var State; wallStart: float; opts: RunOptions) =
  var changed = true
  while changed:
    changed = false
    for i in 0 ..< m.scheduledActions.len:
      if m.scheduledActions[i].active and m.scheduledActions[i].nextTime <= s.t + timeTolerance(m.scheduledActions[i].nextTime, s.t):
        let e = m.scheduledActions[i]
        executeScheduledAction(m, s, e, wallStart, opts)
        if m.scheduledActions[i].repeat:
          if m.scheduledActions[i].remaining > 0:
            dec m.scheduledActions[i].remaining
            if m.scheduledActions[i].remaining == 0:
              m.scheduledActions[i].active = false
          if m.scheduledActions[i].active:
            let kNext = int64(floor((s.t + timeTolerance(s.t, m.scheduledActions[i].startTime) - m.scheduledActions[i].startTime) / m.scheduledActions[i].period)) + 1'i64
            m.scheduledActions[i].nextTime = m.scheduledActions[i].startTime + float(kNext) * m.scheduledActions[i].period
        else:
          m.scheduledActions[i].active = false
        changed = true

proc selectChannel(props: seq[float]; total: float; rng: var Rand): int =
  let r = rng.rand(1.0) * total
  var acc = 0.0
  var lastPositive = -1
  for i, a in props:
    if a <= 0.0:
      continue
    lastPositive = i
    acc += a
    if r < acc:
      return i
  if lastPositive >= 0:
    return lastPositive
  raise newException(ValueError, "selectChannel called without a positive propensity")

proc checkMemoryPolicy(m: Model; s: var State; didSnapshot: var bool): bool =
  ## Parity port from slimmc-copo's checkMemoryPolicy. Returns true if the
  ## run should stop. On limit, "snapshot" triggers one full logical
  ## snapshot (the canonical Storage snapshot columns, all
  ## sharing one snapshot_id via saveSnapshot()) rather than just
  ## the Storage memory snapshot, since the point of a limit-triggered snapshot is to
  ## capture everything available right before the run potentially stops.
  if not m.memoryPolicy.hasLimit: return false
  let mem = estimateMemory(m, s)
  if mem.totalBytes < m.memoryPolicy.limitBytes: return false
  echo "[memory] limit reached: estimated total=", fmtBytes(mem.totalBytes),
       " limit=", fmtBytes(m.memoryPolicy.limitBytes)
  if m.memoryPolicy.snapshotOnLimit and not didSnapshot:
    saveSnapshot(m, s, withChains = true, source = "memory_limit")
    didSnapshot = true
  if m.memoryPolicy.stopOnLimit:
    return true
  result = false

var resultsInterruptRequested {.volatile.}: bool

proc resultsControlCHook() {.noconv.} =
  resultsInterruptRequested = true

proc runSimulation*(m0: Model; opts: RunOptions = RunOptions()) =
  var m = m0
  var s = initState(m)
  var rng = initRand(m.seed)
  let wallStart = epochTime()
  let startedAt = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
  s.startedAt = startedAt
  var stopReason = "t_end"
  var didMemorySnapshot = false
  var runStatus = "completed"

  resultsInterruptRequested = false
  setControlCHook(resultsControlCHook)
  var resultsInitialized = false
  try:
    initOutputFiles(m0, s, opts, startedAt)
    resultsInitialized = true
    if opts.debug:
      initDebugFile(m)
      debugCheckState(m, s, "initial state")

    # Automatic lifecycle messages are deliberately distinct from
    # model-controlled `print_info` progress lines, and the start marker
    # precedes any model action scheduled at t=0.
    printStart(m, opts)

    processDueActions(m, s, wallStart, opts)
    if m.conditionalActions.len > 0:
      checkConditionalActions(m, s, wallStart, opts, "initial")
    if s.stopRequested:
      stopReason = "stop_condition"
      runStatus = "completed"

    while not resultsInterruptRequested and not s.stopRequested and s.t < m.tEnd - timeTolerance(s.t, m.tEnd) and s.kmcEvent < m.maxEvents:
      let tScheduled = nextScheduledActionTime(m.scheduledActions)
      let tLimit = min(tScheduled, m.tEnd)
      let props = computePropensities(m, s)

      var a0 = 0.0
      for i, a in props:
        if a.classify in {fcNan, fcInf, fcNegInf} or a < 0.0:
          raise newException(ValueError, "non-finite or negative propensity at event=" & $s.kmcEvent &
            " channel=" & $i & " value=" & $a)
        a0 += a
      if a0.classify in {fcNan, fcInf, fcNegInf}:
        raise newException(ValueError, "non-finite total propensity at event=" & $s.kmcEvent & " a0=" & $a0)

      if a0 <= 0.0:
        if tLimit < InfTime / 2 and tLimit > s.t + timeTolerance(tLimit, s.t):
          s.t = tLimit
          processDueActions(m, s, wallStart, opts)
          continue
        stopReason = "no_positive_propensity"
        runStatus = "stopped"
        break

      var u = rng.rand(1.0)
      if u <= 0.0:
        u = 1.0e-16
      let tau = -ln(u) / a0
      if tau.classify in {fcNan, fcInf, fcNegInf} or tau <= 0.0:
        raise newException(ValueError, "invalid SSA waiting time at event=" & $s.kmcEvent & " tau=" & $tau & " a0=" & $a0)
      s.sumA0Tau += a0 * tau
      s.sumA0TauSq += (a0 * tau) * (a0 * tau)
      s.countA0Tau += 1

      if s.t + tau > tLimit + timeTolerance(s.t + tau, tLimit):
        s.t = tLimit
        processDueActions(m, s, wallStart, opts)
        continue

      s.t += tau
      let chId = selectChannel(props, a0, rng)
      applyChannel(m, s, chId, rng)
      inc s.stateRevision
      if opts.traceChannelsLimit > 0:
        if s.channelTraceRowsWritten < opts.traceChannelsLimit:
          s.storageTraceKmcEvents.add uint64(s.kmcEvent)
          s.storageTraceTimes.add s.t
          s.storageTraceDt.add tau
          s.storageTraceChannelIds.add uint32(chId)
          s.storageTraceRates.add m.rateValue(m.channels[chId].kId)
          s.storageTracePropensities.add props[chId]
          s.storageTraceTotalPropensities.add a0
          inc s.channelTraceRowsWritten
        else:
          s.channelTraceTruncated = true
      inc s.kmcEvent
      inc s.channelFires[chId]

      if checkMemoryPolicy(m, s, didMemorySnapshot):
        stopReason = "memory_limit"
        runStatus = "stopped"
        break

      if m.conditionalActions.len > 0 and s.kmcEvent mod m.whenCheckEvents == 0:
        checkConditionalActions(m, s, wallStart, opts, "ssa_cadence")
        if s.stopRequested:
          stopReason = "stop_condition"
          runStatus = "completed"
          break

    if s.stopRequested:
      stopReason = "stop_condition"
      runStatus = "completed"

    if not s.stopRequested and s.t >= m.tEnd - timeTolerance(s.t, m.tEnd):
      s.t = m.tEnd
      processDueActions(m, s, wallStart, opts)

    if not s.stopRequested and s.kmcEvent >= m.maxEvents and s.t < m.tEnd - timeTolerance(s.t, m.tEnd):
      stopReason = "max_steps"
      runStatus = "stopped"
      warn(0, "max_steps reached before t_end")

    if resultsInterruptRequested:
      captureStorageV1Snapshot(m, s, "manual", isFinal = false, hasChains = false)
      let wallSeconds = epochTime() - wallStart
      let finishedAt = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
      writeRunInfo(m0, s, opts, "interrupted", "user_interrupt", startedAt, finishedAt, wallSeconds)
      writeRunLogFinal(m, opts, s, "user_interrupt", wallSeconds)
      publishStorageV1(m, s, rsInterrupted, startedAt, finishedAt, wallSeconds, 130, "user_interrupt")
      return

    saveSnapshot(m, s, withChains = true, source = "final", isFinal = true)
    let wallSeconds = epochTime() - wallStart
    let finishedAt = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
    writeRunInfo(m0, s, opts, runStatus, stopReason, startedAt, finishedAt, wallSeconds)
    writeRunLogFinal(m, opts, s, stopReason, wallSeconds)
    if opts.debug:
      debugCheckState(m, s, "final state")
      writeDebugFinal(m, s, stopReason, wallSeconds)
    publishStorageV1(m, s, rsCompleted, startedAt, finishedAt, wallSeconds, 0, stopReason)
    printDone(m, s, opts, runStatus, stopReason, wallSeconds)
  except CatchableError as e:
    if resultsInitialized:
      let wallSeconds = epochTime() - wallStart
      let finishedAt = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
      try:
        captureStorageV1Snapshot(m, s, "manual", isFinal = false, hasChains = false)
        let logFile = open(m.outputDir / "diagnostics" / "run.log", fmAppend)
        logFile.write("Slimmc run failed at " & finishedAt & ": " & e.msg & "\n")
        logFile.close()
        publishStorageV1(m, s, rsFailed, startedAt, finishedAt, wallSeconds, 1, "runtime_error")
      except CatchableError:
        discard
    raise
  finally:
    unsetControlCHook()
