# copo_kmc.nim
# SSA/kMC event engine for slimmc-copo.

import math, random, times, strutils, tables
import copo_types
import copo_stats
import copo_io
import copo_storage
import copo_propensity
export copo_propensity
import ../../common/safe_numeric

proc initState*(m: Model): State =
  result.t = 0.0
  result.kmcEvent = 0
  result.actionNo = 0
  result.stateRevision = 0
  result.parameterStateId = 1
  result.snapshotId = 0
  result.snapshots = @[]
  result.nextChainId = 1
  result.anySnapshotWritten = false
  result.chainsWritten = false
  result.channelTraceRowsWritten = 0
  result.channelTraceTruncated = false
  result.feedMonomerRemainders = newSeq[seq[float]](m.feeds.len)
  result.feedSpeciesRemainders = newSeq[seq[float]](m.feeds.len)
  for fid in 0 ..< m.feeds.len:
    result.feedMonomerRemainders[fid] = newSeq[float](m.monomers.len)
    result.feedSpeciesRemainders[fid] = newSeq[float](m.species.len)
  result.speciesN = newSeq[int64](m.species.len)
  result.monomerN = newSeq[int64](m.monomers.len)
  result.monomerN0 = newSeq[int64](m.monomers.len)
  result.monomerBalance = newSeq[int64](m.monomers.len)
  result.monomerExternalN = newSeq[int64](m.monomers.len)
  result.monomerDosedN = newSeq[int64](m.monomers.len)
  result.speciesExternalN = newSeq[int64](m.species.len)
  result.speciesDosedN = newSeq[int64](m.species.len)
  result.livePools = newSeq[seq[LiveChain]](m.pools.len)
  result.poolEligibleCounts = newSeq[int64](m.pools.len)
  result.poolPropagatableCounts = newSeq[int64](m.pools.len)
  result.poolEligibilityTrackedLen = newSeq[int64](m.pools.len)
  result.channelDepropEligibleCounts = newSeq[int64](m.channels.len)
  result.channelFires = newSeq[int64](m.channels.len)
  result.channelSuccesses = newSeq[int64](m.channels.len)
  result.channelFailures = newSeq[int64](m.channels.len)

  for i, sp in m.species:
    result.speciesN[i] = countFromConc(sp.c0, m.V)
    result.speciesExternalN[i] = result.speciesN[i]
  for i, mo in m.monomers:
    result.monomerN[i] = countFromConc(mo.c0, m.V)
    result.monomerN0[i] = result.monomerN[i]
    result.monomerBalance[i] = result.monomerN[i]
    result.monomerExternalN[i] = result.monomerN[i]


# Eligibility and propensity calculation live in copo_propensity so the
# Storage writer can persist exact snapshot propensities without a cycle.
proc livePoolInvariantErrors*(m: Model; s: State): seq[string] =
  ## Lightweight diagnostics for tests/debugging: live chains must match the
  ## terminal and, where constrained, penultimate metadata of their pool.
  for pi, pool in s.livePools:
    for ci, c in pool:
      if not m.poolAcceptsChain(pi, c):
        result.add "pool " & m.pools[pi].name & " chain_index=" & $ci &
          " last=" & $c.last & " prev=" & $c.prev

proc debugCheckState(m: Model; s: State; opts: RunOptions; source: string) =
  ## Classic-like debug/invariant checkpoint. It writes debug.log when --debug
  ## is enabled and raises on structural invariants that should never be
  ## violated by the SSA engine.
  if not opts.debug: return
  var issues: seq[string] = @[]
  for issue in m.livePoolInvariantErrors(s):
    issues.add "live_pool_invariant: " & issue
  for i, n in s.monomerN:
    if n < 0: issues.add "negative monomer count " & m.monomers[i].name & "=" & $n
  for i, n in s.speciesN:
    if n < 0: issues.add "negative species count " & m.species[i].name & "=" & $n
  let missing = missingEndgroups(m, s)
  let msg = if issues.len == 0:
    "[debug] source=" & source & " event=" & $s.kmcEvent & " t=" & $s.t & " ok missing_endgroups=" & missing.join(",")
  else:
    "[debug] source=" & source & " event=" & $s.kmcEvent & " t=" & $s.t & " issues=" & issues.join(" | ")
  appendDebugLog(m, msg)
  if issues.len > 0:
    raise newException(ValueError, "debugCheckState failed at " & source & ": " & issues.join(" | "))


proc takeAt(pool: var seq[LiveChain]; idx: int): LiveChain =
  result = pool[idx]
  pool[idx] = pool[^1]
  pool.setLen(pool.len - 1)


proc eligibilityCountersCurrent(s: State; poolId: int): bool =
  poolId >= 0 and
  poolId < s.livePools.len and
  poolId < s.poolEligibleCounts.len and
  poolId < s.poolPropagatableCounts.len and
  poolId < s.poolEligibilityTrackedLen.len and
  s.poolEligibilityTrackedLen[poolId] == int64(s.livePools[poolId].len)

proc depropChannelAccepts(m: Model; ch: KmcChannel; c: LiveChain): bool =
  if ch.kind != chMacroDeprop: return false
  if ch.poolOut < 0 or ch.poolOut >= m.poolTerminalMer.len: return false
  let newTerminal = m.poolTerminalMer[ch.poolOut]
  newTerminal >= 0 and c.dp >= 2 and c.last == ch.monomerId and c.prev == newTerminal

proc trackAdd(s: var State; m: Model; poolId: int; c: LiveChain) =
  ## Add through the engine-owned path and update all O(1) eligibility counters.
  let oldLen = s.livePools[poolId].len
  let tracked = s.poolEligibleCounts.len == s.livePools.len and
    s.poolPropagatableCounts.len == s.livePools.len and
    s.poolEligibilityTrackedLen.len == s.livePools.len and
    s.channelDepropEligibleCounts.len == m.channels.len and
    s.poolEligibilityTrackedLen[poolId] == int64(oldLen)
  s.livePools[poolId].add c
  if not tracked:
    # Defensive fallback for manually assembled/mutated State values. Once a
    # pool is untracked it stays scan-backed; applyChannel must not accidentally
    # bless partial counters after removing/adding one of those chains.
    if s.poolEligibilityTrackedLen.len == s.livePools.len:
      s.poolEligibilityTrackedLen[poolId] = -1
    return
  if m.poolAcceptsChain(poolId, c):
    inc s.poolEligibleCounts[poolId]
    if int64(c.dp) < m.dp_max:
      inc s.poolPropagatableCounts[poolId]
    for chId, ch in m.channels:
      if ch.pool1 == poolId and m.depropChannelAccepts(ch, c):
        inc s.channelDepropEligibleCounts[chId]
  s.poolEligibilityTrackedLen[poolId] = int64(s.livePools[poolId].len)

proc trackRemove(s: var State; m: Model; poolId, idx: int): LiveChain =
  ## Remove through the engine-owned path and update all O(1) eligibility counters.
  let c = s.livePools[poolId][idx]
  let tracked = s.eligibilityCountersCurrent(poolId) and
    s.channelDepropEligibleCounts.len == m.channels.len
  result = s.livePools[poolId].takeAt(idx)
  if not tracked:
    if s.poolEligibilityTrackedLen.len == s.livePools.len:
      s.poolEligibilityTrackedLen[poolId] = -1
    return
  if m.poolAcceptsChain(poolId, c):
    dec s.poolEligibleCounts[poolId]
    if int64(c.dp) < m.dp_max:
      dec s.poolPropagatableCounts[poolId]
    for chId, ch in m.channels:
      if ch.pool1 == poolId and m.depropChannelAccepts(ch, c):
        dec s.channelDepropEligibleCounts[chId]
  s.poolEligibilityTrackedLen[poolId] = int64(s.livePools[poolId].len)

proc chooseTwoDistinct(n: int; rng: var Rand): tuple[i: int, j: int] =
  let i = rng.rand(n - 1)
  var j = rng.rand(n - 2)
  if j >= i: inc j
  result = (i: i, j: j)

proc chooseEligibleIndex(m: Model; s: State; poolId: int; rng: var Rand; context: string): int =
  if s.eligibilityCountersCurrent(poolId) and
     s.poolEligibleCounts[poolId] == int64(s.livePools[poolId].len):
    if s.livePools[poolId].len == 0:
      raise newException(ValueError, context & " selected empty live pool")
    return rng.rand(s.livePools[poolId].len - 1)
  let elig = m.eligiblePoolIndices(s, poolId)
  if elig.len == 0:
    raise newException(ValueError, context & " selected pool with no invariant-compatible live chains")
  result = elig[rng.rand(elig.len - 1)]

proc takeTwo(pool: var seq[LiveChain]; i0, j0: int): tuple[a: LiveChain, b: LiveChain] =
  var i = i0
  var j = j0
  if i == j: raise newException(ValueError, "takeTwo got identical indices")
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


proc applySmallRefs(s: var State; refs: seq[SmallRef]; sign: int) =
  ## Applies a list of SmallRef with the given sign (-1 for consuming
  ## reactants, always called; +1 for forming products, called only on
  ## a successful f-coin-flip). Replaces the old applySmallDeltas, which
  ## pre-merged lhs and rhs into one net delta per species -- that
  ## representation couldn't express "reactants always consumed,
  ## products only sometimes formed" (items 79/80/81).
  for r in refs:
    let delta = sign * r.stoich
    case r.kind
    of skSpecies:
      s.speciesN[r.id] += int64(delta)
      if s.speciesN[r.id] < 0:
        raise newException(ValueError, "elementary rxn produced negative species count")
    of skMonomer:
      s.monomerN[r.id] += int64(delta)
      s.monomerBalance[r.id] += int64(delta)
      if s.monomerN[r.id] < 0:
        raise newException(ValueError, "elementary rxn produced negative monomer count")

proc applyChannel*(m: Model; s: var State; chId: int; rng: var Rand) =
  let ch = m.channels[chId]
  let names = m.monomerNames()

  case ch.kind
  of chRxnUni, chRxnBiDiff, chRxnBiSame:
    # Items 79/80/81: reactants are consumed on every firing, regardless
    # of f. Products form only with probability ch.efficiency -- a
    # failed attempt still consumes substrate (a real "unsuccessful
    # reactive collision"), it just produces nothing. Item 82: track
    # successes/failures so fires = successes + failures holds exactly.
    s.applySmallRefs(ch.smallReactants, -1)
    if rng.rand(1.0) <= ch.efficiency:
      s.applySmallRefs(ch.smallProducts, 1)
      s.channelSuccesses[chId] += 1
    else:
      s.channelFailures[chId] += 1

  of chMacroInit:
    s.speciesN[ch.speciesId] -= 1
    s.monomerN[ch.monomerId] -= 1
    let chain = makeLiveChain(m, s.nextChainId, m.species[ch.speciesId].name,
                              ch.monomerId)
    s.nextChainId += 1
    s.trackAdd(m, ch.poolOut, chain)

  of chMacroProp:
    s.monomerN[ch.monomerId] -= 1
    let idx =
      if s.eligibilityCountersCurrent(ch.pool1) and
         s.poolPropagatableCounts[ch.pool1] == int64(s.livePools[ch.pool1].len):
        rng.rand(s.livePools[ch.pool1].len - 1)
      else:
        let propElig = m.eligiblePropIndices(s, ch.pool1)
        if propElig.len == 0:
          raise newException(ValueError, ch.name & " selected pool with no chain below dp_max")
        propElig[rng.rand(propElig.len - 1)]
    var c = s.trackRemove(m, ch.pool1, idx)
    c.pushMonomer(ch.monomerId, m.monomers[ch.monomerId].mw)
    s.trackAdd(m, ch.poolOut, c)

  of chMacroTermC:
    var c1, c2: LiveChain
    if ch.pool1 == ch.pool2:
      var idxA, idxB: int
      if s.eligibilityCountersCurrent(ch.pool1) and
         s.poolEligibleCounts[ch.pool1] == int64(s.livePools[ch.pool1].len):
        if s.livePools[ch.pool1].len < 2:
          raise newException(ValueError, ch.name & " selected pool with fewer than two live chains")
        let ij = chooseTwoDistinct(s.livePools[ch.pool1].len, rng)
        idxA = ij.i; idxB = ij.j
      else:
        let elig = m.eligiblePoolIndices(s, ch.pool1)
        if elig.len < 2:
          raise newException(ValueError, ch.name & " selected pool with fewer than two invariant-compatible live chains")
        let ij = chooseTwoDistinct(elig.len, rng)
        idxA = elig[ij.i]; idxB = elig[ij.j]
      if idxA > idxB:
        c1 = s.trackRemove(m, ch.pool1, idxA)
        c2 = s.trackRemove(m, ch.pool1, idxB)
      else:
        c2 = s.trackRemove(m, ch.pool1, idxB)
        c1 = s.trackRemove(m, ch.pool1, idxA)
    else:
      c1 = s.trackRemove(m, ch.pool1, m.chooseEligibleIndex(s, ch.pool1, rng, ch.name))
      c2 = s.trackRemove(m, ch.pool2, m.chooseEligibleIndex(s, ch.pool2, rng, ch.name))
    s.deadChains.add combineLiveToDead(m, c1, c2, names, m.sequence_mode == "full")

  of chMacroTermD:
    var c1, c2: LiveChain
    if ch.pool1 == ch.pool2:
      var idxA, idxB: int
      if s.eligibilityCountersCurrent(ch.pool1) and
         s.poolEligibleCounts[ch.pool1] == int64(s.livePools[ch.pool1].len):
        if s.livePools[ch.pool1].len < 2:
          raise newException(ValueError, ch.name & " selected pool with fewer than two live chains")
        let ij = chooseTwoDistinct(s.livePools[ch.pool1].len, rng)
        idxA = ij.i; idxB = ij.j
      else:
        let elig = m.eligiblePoolIndices(s, ch.pool1)
        if elig.len < 2:
          raise newException(ValueError, ch.name & " selected pool with fewer than two invariant-compatible live chains")
        let ij = chooseTwoDistinct(elig.len, rng)
        idxA = elig[ij.i]; idxB = elig[ij.j]
      if idxA > idxB:
        c1 = s.trackRemove(m, ch.pool1, idxA)
        c2 = s.trackRemove(m, ch.pool1, idxB)
      else:
        c2 = s.trackRemove(m, ch.pool1, idxB)
        c1 = s.trackRemove(m, ch.pool1, idxA)
    else:
      c1 = s.trackRemove(m, ch.pool1, m.chooseEligibleIndex(s, ch.pool1, rng, ch.name))
      c2 = s.trackRemove(m, ch.pool2, m.chooseEligibleIndex(s, ch.pool2, rng, ch.name))
    if rng.rand(1) == 0:
      s.deadChains.add makeDeadSummary(m, c1, "H", fbTermD_H, names, m.sequence_mode == "full")
      s.deadChains.add makeDeadSummary(m, c2, "U", fbTermD_U, names, m.sequence_mode == "full")
    else:
      s.deadChains.add makeDeadSummary(m, c1, "U", fbTermD_U, names, m.sequence_mode == "full")
      s.deadChains.add makeDeadSummary(m, c2, "H", fbTermD_H, names, m.sequence_mode == "full")

  of chMacroTermX:
    s.speciesN[ch.speciesId] -= 1
    let idx = m.chooseEligibleIndex(s, ch.pool1, rng, ch.name)
    let c = s.trackRemove(m, ch.pool1, idx)
    s.deadChains.add makeDeadSummary(m, c, m.species[ch.speciesId].name, fbTermX, names, m.sequence_mode == "full")

  of chMacroTransfer:
    s.speciesN[ch.speciesId] -= 1
    s.speciesN[ch.speciesOutId] += 1
    let idx = m.chooseEligibleIndex(s, ch.pool1, rng, ch.name)
    let c = s.trackRemove(m, ch.pool1, idx)
    s.deadChains.add makeDeadSummary(m, c, m.species[ch.speciesId].name, fbTransfer, names, m.sequence_mode == "full")

  of chMacroTransferM:
    s.monomerN[ch.monomerId] -= 1
    let idx = m.chooseEligibleIndex(s, ch.pool1, rng, ch.name)
    let c = s.trackRemove(m, ch.pool1, idx)
    s.deadChains.add makeDeadSummary(m, c, "H", fbTransferM, names, m.sequence_mode == "full")
    let trEnd = m.monomers[ch.monomerId].name & "_tr"
    let chain = makeLiveChain(m, s.nextChainId, trEnd, ch.monomerId)
    s.nextChainId += 1
    s.trackAdd(m, ch.pool2, chain)

  of chMacroDeprop:
    let elig = m.eligibleDepropIndices(s, ch)
    let idx = elig[rng.rand(elig.len - 1)]
    var c = s.trackRemove(m, ch.pool1, idx)
    let popped = c.popTerminalMonomer(m.monomers[ch.monomerId].mw)
    if popped != ch.monomerId:
      raise newException(ValueError, "deprop channel selected inconsistent terminal mer")
    s.monomerN[ch.monomerId] += 1
    s.trackAdd(m, ch.poolOut, c)

proc chooseChannel(a: seq[float]; a0: float; rng: var Rand): int =
  let r = rng.rand(1.0) * a0
  var acc = 0.0
  var lastPositive = -1
  for i, x in a:
    if x <= 0.0:
      continue
    lastPositive = i
    acc += x
    if r < acc:
      return i
  if lastPositive >= 0:
    return lastPositive
  raise newException(ValueError, "chooseChannel called without a positive propensity")

proc nextScheduledActionTime(actions: seq[ScheduledAction]): float =
  result = InfTime
  for a in actions:
    if a.active and a.nextTime < result:
      result = a.nextTime

proc actionFail(lineNo: int; msg: string) =
  raise newException(ValueError, "line " & $lineNo & ": " & msg)

proc targetConc(m: Model; s: State; name: string; lineNo: int): tuple[targetKind: string, id: int, value: float] =
  if m.monomerByName.hasKey(name):
    let id = m.monomerByName[name]
    return (targetKind: "monomer", id: id, value: concFromCount(s.monomerN[id], m.V))
  if m.speciesByName.hasKey(name):
    let id = m.speciesByName[name]
    return (targetKind: "species", id: id, value: concFromCount(s.speciesN[id], m.V))
  actionFail(lineNo, "unknown species or monomer: " & name)

proc deltaCountFromConc(c: float; V: float): int64 =
  checkedCountFromConc(c, NA, V, "concentration increment", allowNegative = true)

proc setTargetCount(m: Model; s: var State; targetKind: string; targetId: int; newN: int64) =
  if newN < 0:
    raise newException(ValueError, "target molecule count must be non-negative")
  if targetKind == "monomer":
    let oldN = s.monomerN[targetId]
    let delta = checkedSubInt64(newN, oldN, "set_c/add_c monomer delta")
    let newN0 = checkedAddInt64(s.monomerN0[targetId], delta, "monomer reference count")
    let newBalance = checkedAddInt64(s.monomerBalance[targetId], delta, "monomer balance")
    if newN0 < 0:
      raise newException(ValueError, "set_c/add_c would make monomer reference count negative")
    s.monomerN[targetId] = newN
    s.monomerN0[targetId] = newN0
    s.monomerBalance[targetId] = newBalance
  else:
    s.speciesN[targetId] = newN


proc actionTarget(action: ActionKind; args: seq[string]): string =
  case action
  of eaSetK, eaAddK, eaSetC, eaAddC, eaFeed:
    if args.len >= 1: return args[0]
  of eaSetTemp, eaAddTemp:
    return "temperature"
  else:
    discard
  result = ""

proc actionRequested(args: seq[string]): string =
  result = args.join(" ")

proc actionObservedValue(m: Model; s: State; action: ActionKind; args: seq[string]): string =
  ## Best-effort user-facing value for action trace before/after fields.
  ## Invalid action syntax is still diagnosed by executeAction; this helper is
  ## deliberately non-throwing except for impossible internal state.
  case action
  of eaSetK, eaAddK:
    if args.len >= 1 and m.rateByName.hasKey(args[0]):
      return $m.rateValue(m.rateByName[args[0]])
  of eaSetTemp, eaAddTemp:
    return $m.T
  of eaSetC, eaAddC:
    if args.len >= 1:
      if m.monomerByName.hasKey(args[0]):
        let id = m.monomerByName[args[0]]
        return $concFromCount(s.monomerN[id], m.V)
      if m.speciesByName.hasKey(args[0]):
        let id = m.speciesByName[args[0]]
        return $concFromCount(s.speciesN[id], m.V)
  else:
    discard
  result = ""

proc actionWritesOutput(action: ActionKind): bool =
  case action
  of eaPrint, eaSave, eaSaveChains, eaPrintMemory,
     eaSetK, eaAddK, eaSetTemp, eaAddTemp:
    result = true
  else:
    result = false

proc changesParameterState(action: ActionKind): bool =
  action in {eaSetK, eaAddK, eaSetTemp, eaAddTemp}

proc finalizeStateChange(m: Model; s: var State; action: ActionKind) =
  if not changesParameterState(action):
    return
  s.parameterStateId += 1
  writeParameterState(m, s, s.actionNo)
  captureKineticParameterSetV1(m, s, true, uint64(max(0'i64, s.actionNo - 1)))
  saveSnapshot(m, s, false, "parameter_change")

proc executeAction(m: var Model; s: var State; action: ActionKind; args: seq[string]; lineNo: int; opts: RunOptions; wallStart: float): bool =
  case action
  of eaPrint:
    # Item 57: print "MESSAGE" outputs only the user's message -- no
    # "[action] t=..." prefix (previously present here, unlike homo's
    # printMarker(), which has always written just the bare message).
    stdout.writeLine(args[0])
    result = false
  of eaPrintInfo:
    printStatus(m, s, wallStart, opts)
    result = false
  of eaSave:
    if args.len != 0: actionFail(lineNo, "save takes no arguments")
    saveSnapshot(m, s, false, "save")
    result = false
  of eaSaveChains:
    if args.len != 0: actionFail(lineNo, "save_chains takes no arguments")
    saveSnapshot(m, s, true, "save_chains")
    result = false
  of eaStop:
    if args.len != 0: actionFail(lineNo, "stop takes no arguments")
    result = false
  of eaPrintMemory:
    if args.len != 0: actionFail(lineNo, "print_memory takes no arguments")
    printMemory(m, s)
    result = false
  of eaSetK:
    if args.len != 2: actionFail(lineNo, "set_k syntax: set_k rate value")
    if not m.rateByName.hasKey(args[0]): actionFail(lineNo, "unknown rate: " & args[0])
    let kId = m.rateByName[args[0]]
    let newK = parseFloat(args[1])
    if newK < 0.0: actionFail(lineNo, "set_k requires rate value >= 0")
    m.rates[kId].kind = rkFixed
    m.rates[kId].kConst = newK
    result = true
  of eaAddK:
    if args.len != 2: actionFail(lineNo, "add_k syntax: add_k rate increment")
    if not m.rateByName.hasKey(args[0]): actionFail(lineNo, "unknown rate: " & args[0])
    let kId = m.rateByName[args[0]]
    let newK = m.rateValue(kId) + parseFloat(args[1])
    if newK < 0.0: actionFail(lineNo, "add_k would make negative rate constant")
    m.rates[kId].kind = rkFixed
    m.rates[kId].kConst = newK
    result = true
  of eaSetTemp:
    if args.len != 1: actionFail(lineNo, "set_temp syntax: set_temp value")
    let newT = parseFloat(args[0])
    if newT <= 0.0: actionFail(lineNo, "set_temp requires temperature > 0")
    m.T = newT
    result = true
  of eaAddTemp:
    if args.len != 1: actionFail(lineNo, "add_temp syntax: add_temp increment")
    let newT = m.T + parseFloat(args[0])
    if newT <= 0.0: actionFail(lineNo, "add_temp would make non-positive temperature")
    m.T = newT
    result = true
  of eaSetC:
    if args.len != 2: actionFail(lineNo, "set_c syntax: set_c species_or_monomer value")
    let target = targetConc(m, s, args[0], lineNo)
    let warning = "WARNING: line " & $lineNo & ": set_c forces the concentration of '" & args[0] & "'. Its physical material balance is invalid after this action."
    stdout.writeLine(warning)
    appendRunLog(m, warning)
    let c = parseFloat(args[1])
    if c < 0.0: actionFail(lineNo, "set_c requires concentration >= 0")
    setTargetCount(m, s, target.targetKind, target.id, countFromConc(c, m.V))
    result = true
  of eaAddC:
    if args.len != 2: actionFail(lineNo, "add_c syntax: add_c species_or_monomer increment")
    let target = targetConc(m, s, args[0], lineNo)
    let dc = parseFloat(args[1])
    let dn = deltaCountFromConc(dc, m.V)
    let oldN = if target.targetKind == "monomer": s.monomerN[target.id] else: s.speciesN[target.id]
    let newN = checkedAddInt64(oldN, dn, "add_c molecule count")
    if newN < 0: actionFail(lineNo, "add_c would make negative molecule count")
    setTargetCount(m, s, target.targetKind, target.id, newN)
    if target.targetKind == "monomer":
      s.monomerExternalN[target.id] = checkedAddInt64(s.monomerExternalN[target.id], dn, "add_c external monomer balance")
    else:
      s.speciesExternalN[target.id] = checkedAddInt64(s.speciesExternalN[target.id], dn, "add_c external species balance")
    result = true

  of eaFeed:
    if args.len notin [2, 3]: actionFail(lineNo, "feed action syntax: feed NAME VOLUME [L|l|mL|ml|ML]; without a unit VOLUME is in L")
    if not m.feedByName.hasKey(args[0]): actionFail(lineNo, "unknown feed: " & args[0])
    let doseValue = parseFloat(args[1])
    var doseMl: float
    if args.len == 2:
      doseMl = doseValue * 1000.0
    else:
      let unit = args[2].toLowerAscii()
      case unit
      of "l": doseMl = doseValue * 1000.0
      of "ml": doseMl = doseValue
      else: actionFail(lineNo, "feed volume unit must be L or mL (case-insensitive)")
    if doseMl <= 0.0: actionFail(lineNo, "feed volume must be > 0")
    if m.currentVolumeMl <= 0.0: actionFail(lineNo, "feed requires param init_volume > 0")
    let fid = m.feedByName[args[0]]
    let oldV = m.V
    let newPhysical = m.currentVolumeMl + doseMl
    let newV = oldV * newPhysical / m.currentVolumeMl
    let deltaV = newV - oldV
    for mid in 0 ..< m.monomers.len:
      let exact = m.feeds[fid].monomerConcentrations[mid] * NA * deltaV + s.feedMonomerRemainders[fid][mid]
      if exact < 0.0 or exact > float(high(int64)): actionFail(lineNo, "feed molecule count overflow for " & m.monomers[mid].name)
      let dn = int64(floor(exact))
      s.feedMonomerRemainders[fid][mid] = exact - float(dn)
      let newN = checkedAddInt64(s.monomerN[mid], dn, "feed monomer count")
      setTargetCount(m, s, "monomer", mid, newN)
      s.monomerExternalN[mid] = checkedAddInt64(s.monomerExternalN[mid], dn, "feed external monomer balance")
      s.monomerDosedN[mid] = checkedAddInt64(s.monomerDosedN[mid], dn, "feed dosed monomer balance")
    for sid in 0 ..< m.species.len:
      let feedC = (if sid < m.feeds[fid].speciesConcentrations.len: m.feeds[fid].speciesConcentrations[sid] else: 0.0)
      let exact = feedC * NA * deltaV + s.feedSpeciesRemainders[fid][sid]
      if exact < 0.0 or exact > float(high(int64)): actionFail(lineNo, "feed molecule count overflow for " & m.species[sid].name)
      let dn = int64(floor(exact))
      s.feedSpeciesRemainders[fid][sid] = exact - float(dn)
      s.speciesN[sid] = checkedAddInt64(s.speciesN[sid], dn, "feed species count")
      s.speciesExternalN[sid] = checkedAddInt64(s.speciesExternalN[sid], dn, "feed external species balance")
      s.speciesDosedN[sid] = checkedAddInt64(s.speciesDosedN[sid], dn, "feed dosed species balance")
    m.currentVolumeMl = newPhysical
    m.V = newV
    result = true
  if true:
    appendRunLog(m, "[action] event=" & $s.kmcEvent & " t=" & $s.t & " action=" & actionKindName(action) & " args=" & args.join(" "))

proc totalConversion(s: State): float =
  var cur = 0'i64
  var exp = 0'i64
  for x in s.monomerN: cur += x
  for x in s.monomerN0: exp += x
  if exp <= 0: return 0.0
  result = 1.0 - float(cur) / float(exp)

proc conditionalValue(m: Model; s: State; a: AtomicCondition): float =
  case a.observable
  of woTotalConversion:
    result = totalConversion(s)
  of woMonomerConversion:
    if s.monomerN0[a.targetId] <= 0: result = 0.0
    else: result = 1.0 - float(s.monomerN[a.targetId]) / float(s.monomerN0[a.targetId])
  of woSpeciesConc:
    result = concFromCount(s.speciesN[a.targetId], m.V)
  of woMonomerConc:
    result = concFromCount(s.monomerN[a.targetId], m.V)

proc conditionalIsTrue(value: float; a: AtomicCondition): bool =
  case a.comparison
  of coGreater: result = value > a.threshold
  of coLess: result = value < a.threshold

proc conditionText(m: Model; e: ConditionalAction): string =
  var parts: seq[string] = @[]
  for a in e.conditions:
    parts.add conditionalObservableName(m, a) & " " & comparisonName(a.comparison) & " " & $a.threshold
  result = parts.join(" and ")

proc conditionValues(values: seq[float]): string =
  var parts: seq[string] = @[]
  for v in values: parts.add $v
  result = parts.join(" and ")

proc checkConditionalActions(m: var Model; s: var State; opts: RunOptions; source: string; wallStart: float) =
  var fired = true
  while fired:
    fired = false
    for i in 0 ..< m.conditionalActions.len:
      if not m.conditionalActions[i].active: continue
      let e = m.conditionalActions[i]
      var values: seq[float] = @[]
      var matched = true
      for a in e.conditions:
        let value = conditionalValue(m, s, a)
        values.add value
        if not conditionalIsTrue(value, a): matched = false
      if not matched: continue
      m.conditionalActions[i].active = false
      s.actionNo += 1
      let beforeValue = actionObservedValue(m, s, e.action, e.args)
      let stateChanged = executeAction(m, s, e.action, e.args, e.lineNo, opts, wallStart)
      if stateChanged:
        s.stateRevision += 1
        finalizeStateChange(m, s, e.action)
      if e.action == eaStop:
        s.stopRequested = true
        s.stopLineNo = e.lineNo
        s.stopCheckSource = source
        s.stopConditions = e.conditions
        s.stopActualValues = values
      let afterValue = actionObservedValue(m, s, e.action, e.args)
      writeActionTrace(m, s, "when", "", source, conditionText(m, e),
        (if e.conditions.len > 1: "and" else: comparisonName(e.conditions[0].comparison)),
        (if e.conditions.len > 1: conditionText(m, e) else: $e.conditions[0].threshold),
        conditionValues(values), e.action,
        actionTarget(e.action, e.args), actionRequested(e.args), beforeValue, afterValue, "ok",
        e.lineNo, stateChanged, actionWritesOutput(e.action))
      captureActionV1(m, s, e.lineNo, "when", "", e.conditions, values, e.action, e.args,
        beforeValue, afterValue, stateChanged, actionWritesOutput(e.action),
        (if e.action == eaPrint and e.args.len > 0: e.args[0] else: ""))
      if true:
        appendRunLog(m, "[when] line=" & $e.lineNo & " source=" & source & " condition=" & conditionText(m, e) & " values=" & conditionValues(values) & " action=" & actionKindName(e.action))
      fired = true
      discard stateChanged

proc processDueActions(m: var Model; s: var State; opts: RunOptions; wallStart: float): bool =
  var changed = true
  while changed:
    changed = false
    for i in 0 ..< m.scheduledActions.len:
      if m.scheduledActions[i].active and m.scheduledActions[i].nextTime <= s.t + timeTolerance(m.scheduledActions[i].nextTime, s.t):
        let e = m.scheduledActions[i]
        s.actionNo += 1
        let beforeValue = actionObservedValue(m, s, e.action, e.args)
        let stateChanged = executeAction(m, s, e.action, e.args, e.lineNo, opts, wallStart)
        if stateChanged:
          s.stateRevision += 1
          finalizeStateChange(m, s, e.action)
        let afterValue = actionObservedValue(m, s, e.action, e.args)
        result = true
        writeActionTrace(m, s, (if e.repeat: "every" else: "at"), $e.nextTime,
          "scheduled_action", "", "", "", "", e.action,
          actionTarget(e.action, e.args), actionRequested(e.args), beforeValue, afterValue, "ok",
          e.lineNo, stateChanged, actionWritesOutput(e.action))
        captureActionV1(m, s, e.lineNo, (if e.repeat: "every" else: "at"), $e.nextTime,
          @[], @[], e.action, e.args, beforeValue, afterValue, stateChanged, actionWritesOutput(e.action),
          (if e.action == eaPrint and e.args.len > 0: e.args[0] else: ""))
        if true:
          appendRunLog(m, "[scheduled] line=" & $e.lineNo & " kind=" & (if e.repeat: "every" else: "at") & " time=" & $e.nextTime & " action=" & actionKindName(e.action))
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
        checkConditionalActions(m, s, opts, "scheduled_action", wallStart)
        changed = true

proc checkMemoryPolicy(m: Model; s: var State; didSnapshot: var bool): bool =
  if not m.memoryPolicy.hasLimit: return false
  let mem = estimateMemory(m, s)
  if mem.totalBytes < m.memoryPolicy.limitBytes: return false
  echo "[memory] limit reached: estimated total=", fmtBytes(mem.totalBytes),
       " limit=", fmtBytes(m.memoryPolicy.limitBytes)
  if m.memoryPolicy.snapshotOnLimit and not didSnapshot:
    saveSnapshot(m, s, true, "memory_limit")
    didSnapshot = true
  if m.memoryPolicy.stopOnLimit:
    return true
  result = false

var resultsInterruptRequested {.volatile.}: bool

proc resultsControlCHook() {.noconv.} =
  resultsInterruptRequested = true

proc runSimulation*(m: Model; opts: RunOptions = RunOptions()) =
  var model = m
  var s = initState(model)
  var rng = initRand(model.seed)
  var wallStart = epochTime()
  var didMemorySnapshot = false
  var runStatus = "running"
  var stopReason = ""

  resultsInterruptRequested = false
  setControlCHook(resultsControlCHook)

  model.resetOutputDir()
  let startedAt = writeStartInfo(model, opts)
  model.startedAt = startedAt
  initStorageV1(model, s)
  initDiagnostics(model, opts)
  initParameterStates(model, s)
  debugCheckState(model, s, opts, "initial")
  if true:
    appendRunLog(model, "[start] run_id=" & model.modelStem & " t_end=" & $model.t_end & " max_steps=" & $model.max_steps)

  discard processDueActions(model, s, opts, wallStart)
  checkConditionalActions(model, s, opts, "initial", wallStart)
  if s.stopRequested:
    runStatus = "completed"
    stopReason = "stop_condition"
  printStatus(model, s, wallStart, opts)
  # NOTE (family parity): classic slimmc never
  # writes state/moments/etc. unless the model explicitly asks for it via a
  # `save`/`save_chains` action. The unconditional snapshot that used to be
  # taken here caused a real duplicate-row bug: a model with `at 0 save` and
  # `t_end 0.0` produced duplicate t=0 snapshot rows
  # (one from the explicit action, one from this call, one from the mirror
  # call at the end of the run). Removed so that Storage snapshot tables
  # only ever gain a row when the model author writes `save`/`save_chains`
  # themselves -- exactly like classic slimmc. If you want a guaranteed
  # initial snapshot, write `at 0 save` (and `at 0 save_chains` if needed)
  # in the .model file.

  while not resultsInterruptRequested and not s.stopRequested and s.t < model.t_end - timeTolerance(s.t, model.t_end) and s.kmcEvent < model.max_steps:
    discard processDueActions(model, s, opts, wallStart)
    let tScheduled = nextScheduledActionTime(model.scheduledActions)
    let tLimit = min(tScheduled, model.t_end)
    let a = computePropensities(model, s)
    var a0 = 0.0
    for i, x in a:
      if x.classify in {fcNan, fcInf, fcNegInf} or x < 0.0:
        raise newException(ValueError, "non-finite or negative propensity at event=" & $s.kmcEvent &
          " channel=" & $i & " value=" & $x)
      a0 += x
    if a0.classify in {fcNan, fcInf, fcNegInf}:
      raise newException(ValueError, "non-finite total propensity at event=" & $s.kmcEvent & " a0=" & $a0)

    if a0 <= 0.0:
      if tLimit < InfTime / 2 and tLimit > s.t + timeTolerance(tLimit, s.t):
        s.t = tLimit
        discard processDueActions(model, s, opts, wallStart)
        checkConditionalActions(model, s, opts, "time_barrier_no_channels", wallStart)
        continue
      runStatus = "stopped"
      stopReason = "no_active_channels"
      echo "[run] stopped normally: no active channels at event=", s.kmcEvent, " t=", s.t
      if true:
        appendRunLog(model, "[stop] no active channels event=" & $s.kmcEvent & " t=" & $s.t)
      break

    let r1 = max(rng.rand(1.0), 1.0e-16)
    let tau = -ln(r1) / a0
    if tau.classify in {fcNan, fcInf, fcNegInf} or tau <= 0.0:
      raise newException(ValueError, "invalid SSA waiting time at event=" & $s.kmcEvent & " tau=" & $tau & " a0=" & $a0)
    s.sumA0Tau += a0 * tau
    s.sumA0TauSq += (a0 * tau) * (a0 * tau)
    s.countA0Tau += 1
    if s.t + tau > tLimit + timeTolerance(s.t + tau, tLimit):
      s.t = tLimit
      discard processDueActions(model, s, opts, wallStart)
      checkConditionalActions(model, s, opts, "time_barrier", wallStart)
      if s.stopRequested:
        runStatus = "completed"
        stopReason = "stop_condition"
        break
      if s.t >= model.t_end - timeTolerance(s.t, model.t_end):
        break
      continue

    s.t += tau
    let chId = chooseChannel(a, a0, rng)
    applyChannel(model, s, chId, rng)
    s.stateRevision += 1
    s.channelFires[chId] += 1
    s.kmcEvent += 1
    debugCheckState(model, s, opts, "after_channel:" & model.channels[chId].name)
    if opts.traceChannelsLimit > 0:
      if s.channelTraceRowsWritten < opts.traceChannelsLimit:
        s.storageTraceKmcEvents.add uint64(s.kmcEvent)
        s.storageTraceTimes.add s.t
        s.storageTraceDt.add tau
        s.storageTraceChannelIds.add uint32(chId)
        s.storageTraceRates.add model.rateValue(model.channels[chId].kId)
        s.storageTracePropensities.add a[chId]
        s.storageTraceTotalPropensities.add a0
        inc s.channelTraceRowsWritten
      else:
        s.channelTraceTruncated = true

    if model.conditionalActions.len > 0 and s.kmcEvent mod model.whenCheckEvents == 0:
      checkConditionalActions(model, s, opts, "ssa_cadence", wallStart)
      if s.stopRequested:
        runStatus = "completed"
        stopReason = "stop_condition"
        break

    if opts.debug and (s.kmcEvent mod 10000 == 0):
      printMemory(model, s)

    if checkMemoryPolicy(model, s, didMemorySnapshot):
      runStatus = "stopped"
      stopReason = "memory_limit"
      break

  if s.stopRequested:
    runStatus = "completed"
    stopReason = "stop_condition"

  if not s.stopRequested and s.t >= model.t_end - timeTolerance(s.t, model.t_end):
    s.t = model.t_end
    if runStatus == "running":
      runStatus = "completed"
      stopReason = "t_end"
    discard processDueActions(model, s, opts, wallStart)
    checkConditionalActions(model, s, opts, "final", wallStart)
    if s.stopRequested:
      runStatus = "completed"
      stopReason = "stop_condition"
  elif not s.stopRequested and s.kmcEvent >= model.max_steps and runStatus == "running":
    runStatus = "stopped"
    stopReason = "max_steps"

  if runStatus == "running":
    runStatus = "stopped"
    stopReason = "loop_exit"

  if resultsInterruptRequested:
    captureStorageV1Snapshot(model, s, "manual", false, false, computePropensities(model, s))
    let wallSeconds = epochTime() - wallStart
    appendRunLog(model, "[done] event=" & $s.kmcEvent & " t=" & $s.t & " status=interrupted reason=user_interrupt output=" & model.output_dir)
    finalizeInterruptedStorageV1(model, s, wallSeconds, 130)
    writeRunInfo(model, s, opts, "interrupted", "user_interrupt", wallSeconds, startedAt)
    echo "[interrupted] event=", s.kmcEvent, " t=", s.t, " output=", model.output_dir
    unsetControlCHook()
    return

  let finalPrintedByAction = processDueActions(model, s, opts, wallStart)
  if not finalPrintedByAction:
    printStatus(model, s, wallStart, opts)
  debugCheckState(model, s, opts, "final")
  saveSnapshot(model, s, true, "final", true)
  let wallSeconds = epochTime() - wallStart
  appendRunLog(model, "[done] event=" & $s.kmcEvent & " t=" & $s.t & " status=" & runStatus & " reason=" & stopReason & " output=" & model.output_dir)
  finalizeStorageV1(model, s, wallSeconds, stopReason)
  writeRunInfo(model, s, opts, runStatus, stopReason, wallSeconds, startedAt)
  echo "[done] event=", s.kmcEvent, " t=", s.t, " output=", model.output_dir
  unsetControlCHook()
