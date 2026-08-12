import unittest, random, math, tables, algorithm, sequtils
import slimmc_types
import slimmc_kmc

proc phaseAModel(): Model =
  result = initModel()
  result.V = 1.0e-18
  result.T = 298.15
  result.dpMax = 100
  result.species = @[
    SpeciesDef(name: "M", kind: skMonomer, c0: 0.0, mw: 100.0, hasMw: true)
  ]
  result.speciesByName["M"] = 0
  result.monomerId = 0
  result.pools = @[
    PoolDef(name: "P", kind: pkActive),
    PoolDef(name: "Q", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolByName["P"] = 0
  result.poolByName["Q"] = 1
  result.poolByName["D"] = 2
  result.rates = @[
    RateDef(name: "k", kind: rkFixed, kConst: 2.0)
  ]
  result.rateByName["k"] = 0
  discard result.ensureBuiltinEg("R", 68.0, "test")
  discard result.ensureBuiltinEg("S", 44.0, "test")
  discard result.ensureBuiltinEg("M_tr", -1.008, "test")

proc stateFor(m: Model; freeM: int64; p: seq[Chain]; q: seq[Chain] = @[]; dead: seq[Chain] = @[]): State =
  result = initState(m)
  result.n[0] = freeM
  result.pools[0] = p
  result.pools[1] = q
  result.pools[2] = dead
  result.mExpected = calcMTotal(m, result)
  result.mBalance = result.mExpected

proc liveChain(m: Model; dp: int; left = "R"; formedBy = fbInit): Chain =
  result = Chain(eg1: m.egByName[left], dp: dp, eg2: m.egActive, formedBy: formedBy)

proc chainCount(s: State; kind: PoolKind): int =
  for pid, pool in s.pools:
    if pid < s.pools.len and pool.len >= 0:
      discard
  # Model is fixed in these tests: P,Q active; D dead.
  if kind == pkActive: result = s.pools[0].len + s.pools[1].len
  else: result = s.pools[2].len

proc assertCoreInvariants(m: Model; s: State; expectedTotal: int64) =
  check calcMTotal(m, s) == expectedTotal
  for n in s.n:
    check n >= 0
  for pid, pool in s.pools:
    for ch in pool:
      check ch.dp >= 1
      if m.pools[pid].kind == pkActive:
        check ch.eg2 == m.egActive
      else:
        check ch.eg2 != m.egActive

proc closeEnough(a, b: float; rtol = 1.0e-12): bool =
  if a == b: return true
  abs(a - b) <= max(abs(a), abs(b)) * rtol

suite "homo chemistry phase A detailed":
  test "H01 propagation exact propensity, bookkeeping, and dp_max exclusion":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0)]
    var s = stateFor(m, 10, @[liveChain(m, 3), liveChain(m, 100)])
    let total0 = calcMTotal(m, s)
    let expected = 2.0 / (NA * m.V) * 1.0 * 10.0
    check closeEnough(computePropensities(m, s)[0], expected)
    var rng = initRand(101)
    applyChannel(m, s, 0, rng)
    check s.n[0] == 9
    check s.pools[0][0].dp == 4
    check s.pools[0][1].dp == 100
    check chainCount(s, pkActive) == 2
    check chainCount(s, pkDead) == 0
    assertCoreInvariants(m, s, total0)

  test "H01 propensity scales with k, monomer, eligible chains, and inverse volume":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0)]
    let s1 = stateFor(m, 5, @[liveChain(m, 2)])
    let a1 = computePropensities(m, s1)[0]
    var s2 = stateFor(m, 10, @[liveChain(m, 2), liveChain(m, 3)])
    let a2 = computePropensities(m, s2)[0]
    check closeEnough(a2 / a1, 4.0)
    m.rates[0].kConst = 6.0
    check closeEnough(computePropensities(m, s2)[0] / a2, 3.0)
    m.V *= 2.0
    check closeEnough(computePropensities(m, s2)[0] / a2, 1.5)

  test "H01 zero propensity for no monomer or all chains at dp_max":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0)]
    check computePropensities(m, stateFor(m, 0, @[liveChain(m, 2)]))[0] == 0.0
    check computePropensities(m, stateFor(m, 10, @[liveChain(m, 100), liveChain(m, 100)]))[0] == 0.0

  test "H02 deprop exact propensity, monomer return, and DP-one barrier":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0)]
    var s = stateFor(m, 0, @[liveChain(m, 1), liveChain(m, 4), liveChain(m, 7)])
    let total0 = calcMTotal(m, s)
    check computePropensities(m, s)[0] == 4.0
    var rng = initRand(102)
    applyChannel(m, s, 0, rng)
    check s.n[0] == 1
    check s.pools[0].mapIt(it.dp).foldl(a + b) == 11
    check chainCount(s, pkActive) == 3
    assertCoreInvariants(m, s, total0)
    check computePropensities(m, stateFor(m, 0, @[liveChain(m, 1)]))[0] == 0.0

  test "H02 propensity scales exactly with k and eligible-chain count only":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0)]
    let a1 = computePropensities(m, stateFor(m, 999, @[liveChain(m, 1), liveChain(m, 2)]))[0]
    let a2 = computePropensities(m, stateFor(m, 0, @[liveChain(m, 2), liveChain(m, 3), liveChain(m, 1)]))[0]
    check a1 == 2.0
    check a2 == 4.0
    m.rates[0].kConst = 10.0
    check computePropensities(m, stateFor(m, 0, @[liveChain(m, 2), liveChain(m, 3)]))[0] == 20.0

  test "H02 prop then deprop is an exact reversible bookkeeping pair":
    var m = phaseAModel()
    m.channels = @[
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0)
    ]
    var s = stateFor(m, 3, @[liveChain(m, 5)])
    let initial = s
    let total0 = calcMTotal(m, s)
    var rng = initRand(106)
    applyChannel(m, s, 0, rng)
    applyChannel(m, s, 1, rng)
    check s.n == initial.n
    check s.pools == initial.pools
    assertCoreInvariants(m, s, total0)

  test "H03 same-pool combination exact propensity and conservation":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "term_c", kind: chMacroTermC, kId: 0, pool1: 0, pool2: 0, poolOut: 2)]
    var s = stateFor(m, 7, @[liveChain(m, 3), liveChain(m, 8), liveChain(m, 12)])
    let total0 = calcMTotal(m, s)
    let expected = 2.0 / (NA * m.V) * 3.0 * 2.0
    check closeEnough(computePropensities(m, s)[0], expected)
    var rng = initRand(103)
    applyChannel(m, s, 0, rng)
    check chainCount(s, pkActive) == 1
    check chainCount(s, pkDead) == 1
    check s.pools[2][0].formedBy == fbTermC
    check s.pools[2][0].eg2 != m.egActive
    check s.pools[0][0].dp + s.pools[2][0].dp == 23
    assertCoreInvariants(m, s, total0)

  test "H03 cross-pool combination uses n1*n2 and preserves both left ends":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "term_c", kind: chMacroTermC, kId: 0, pool1: 0, pool2: 1, poolOut: 2)]
    var s = stateFor(m, 0, @[liveChain(m, 4, "R")], @[liveChain(m, 9, "S")])
    let total0 = calcMTotal(m, s)
    check closeEnough(computePropensities(m, s)[0], 2.0 / (NA * m.V))
    var rng = initRand(203)
    applyChannel(m, s, 0, rng)
    check chainCount(s, pkActive) == 0
    check s.pools[2].len == 1
    check s.pools[2][0].dp == 13
    check s.pools[2][0].eg1 == m.egByName["R"]
    check s.pools[2][0].eg2 == m.egByName["S"]
    assertCoreInvariants(m, s, total0)

  test "H03 requires two distinct chains in a same-pool reaction":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "term_c", kind: chMacroTermC, kId: 0, pool1: 0, pool2: 0, poolOut: 2)]
    check computePropensities(m, stateFor(m, 0, @[liveChain(m, 5)]))[0] == 0.0

  test "H04 disproportionation preserves individual DPs and creates exactly H plus U":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "term_d", kind: chMacroTermD, kId: 0, pool1: 0, pool2: 0, poolOut: 2)]
    var s = stateFor(m, 5, @[liveChain(m, 4), liveChain(m, 9)])
    let total0 = calcMTotal(m, s)
    check closeEnough(computePropensities(m, s)[0], 2.0 / (NA * m.V) * 2.0)
    var rng = initRand(104)
    applyChannel(m, s, 0, rng)
    check chainCount(s, pkActive) == 0
    check chainCount(s, pkDead) == 2
    var dps = @[s.pools[2][0].dp, s.pools[2][1].dp]
    dps.sort()
    check dps == @[4, 9]
    check s.pools[2].countIt(it.formedBy == fbTermD_H and it.eg2 == m.egH) == 1
    check s.pools[2].countIt(it.formedBy == fbTermD_U and it.eg2 == m.egU) == 1
    assertCoreInvariants(m, s, total0)

  test "H04 H/U assignment is statistically symmetric without changing chemistry":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "term_d", kind: chMacroTermD, kId: 0, pool1: 0, pool2: 1, poolOut: 2)]
    var hOnR = 0
    const nRep = 4000
    for seed in 0 ..< nRep:
      var s = stateFor(m, 0, @[liveChain(m, 4, "R")], @[liveChain(m, 9, "S")])
      var rng = initRand(seed + 1000)
      applyChannel(m, s, 0, rng)
      for ch in s.pools[2]:
        if ch.eg1 == m.egByName["R"] and ch.formedBy == fbTermD_H:
          inc hOnR
    let fraction = float(hOnR) / float(nRep)
    check abs(fraction - 0.5) < 0.04

  test "H05 transfer-to-monomer exact propensity and complete bookkeeping":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "transfer_m", kind: chMacroTransferM, kId: 0, pool1: 0, pool2: 0, poolOut: 2, sp1: 0)]
    var s = stateFor(m, 6, @[liveChain(m, 7)])
    let total0 = calcMTotal(m, s)
    check closeEnough(computePropensities(m, s)[0], 2.0 / (NA * m.V) * 6.0)
    var rng = initRand(105)
    applyChannel(m, s, 0, rng)
    check s.n[0] == 5
    check chainCount(s, pkActive) == 1
    check chainCount(s, pkDead) == 1
    check s.pools[0][0].dp == 1
    check s.pools[0][0].formedBy == fbTransferM
    check s.pools[0][0].eg1 == m.egByName["M_tr"]
    check s.pools[2][0].dp == 7
    check s.pools[2][0].formedBy == fbTransferM
    check s.pools[2][0].eg2 == m.egH
    assertCoreInvariants(m, s, total0)

  test "H05 transfer propensity scales with live chains, monomer, k, and inverse volume":
    var m = phaseAModel()
    m.channels = @[KmcChannel(name: "transfer_m", kind: chMacroTransferM, kId: 0, pool1: 0, pool2: 0, poolOut: 2, sp1: 0)]
    let a1 = computePropensities(m, stateFor(m, 3, @[liveChain(m, 2)]))[0]
    let a2 = computePropensities(m, stateFor(m, 6, @[liveChain(m, 2), liveChain(m, 3)]))[0]
    check closeEnough(a2 / a1, 4.0)
    m.rates[0].kConst = 4.0
    check closeEnough(computePropensities(m, stateFor(m, 6, @[liveChain(m, 2), liveChain(m, 3)]))[0] / a2, 2.0)
    m.V *= 4.0
    check closeEnough(computePropensities(m, stateFor(m, 6, @[liveChain(m, 2), liveChain(m, 3)]))[0] / a2, 0.5)

  test "H06 invariant helper detects valid states after every phase-A channel":
    var m = phaseAModel()
    let cases = @[
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0),
      KmcChannel(name: "term_c", kind: chMacroTermC, kId: 0, pool1: 0, pool2: 0, poolOut: 2),
      KmcChannel(name: "term_d", kind: chMacroTermD, kId: 0, pool1: 0, pool2: 0, poolOut: 2),
      KmcChannel(name: "transfer_m", kind: chMacroTransferM, kId: 0, pool1: 0, pool2: 0, poolOut: 2, sp1: 0)
    ]
    for i, ch in cases:
      m.channels = @[ch]
      var s: State
      case ch.kind
      of chMacroProp: s = stateFor(m, 5, @[liveChain(m, 3)])
      of chMacroDeprop: s = stateFor(m, 0, @[liveChain(m, 3)])
      of chMacroTermC, chMacroTermD: s = stateFor(m, 0, @[liveChain(m, 3), liveChain(m, 4)])
      of chMacroTransferM: s = stateFor(m, 5, @[liveChain(m, 3)])
      else: discard
      let total0 = calcMTotal(m, s)
      var rng = initRand(600 + i)
      applyChannel(m, s, 0, rng)
      assertCoreInvariants(m, s, total0)
