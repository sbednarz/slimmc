import unittest, random
import copo_types
import copo_kmc
import copo_sequence

proc modelForEngineHardening063(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 4
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.pools = @[
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "PB", kind: pkActive),
    PoolDef(name: "PAB", kind: pkActive),
    PoolDef(name: "PBA", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  # Terminal pools intentionally use mixed penultimate metadata.  They must
  # accept AA* and BA* in PA, or AB* and BB* in PB, while still rejecting
  # chains whose terminal mer does not match the pool.
  result.poolTerminalMer = @[0, 1, 1, 0, -1]
  result.poolPenultimateMer = @[-2, -2, 0, 1, -1]
  result.rates = @[
    RateDef(name: "kp_pa_b", kConst: 1.0),
    RateDef(name: "ktc_ab_ba", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "prop_PA_B", kind: chMacroProp, kId: 0,
               pool1: 0, monomerId: 1, poolOut: 2),
    KmcChannel(name: "term_c_PAB_PBA", kind: chMacroTermC, kId: 1,
               pool1: 2, pool2: 3, poolOut: 4)
  ]

proc v063LiveAA(id: int64 = 1): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  result.pushMonomer(0, 100.0)

proc v063LiveBA(id: int64 = 2): LiveChain =
  result = makeLiveChain(id, "R", 1, 128.0, 2)
  result.pushMonomer(0, 100.0)

proc v063LiveAB(id: int64 = 3): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  result.pushMonomer(1, 128.0)

proc v063LiveBB(id: int64 = 4): LiveChain =
  result = makeLiveChain(id, "R", 1, 128.0, 2)
  result.pushMonomer(1, 128.0)

suite "engine regression hardening v0.6.3":
  test "terminal pools with mixed penultimate accept all matching-terminal chains":
    var m = modelForEngineHardening063()
    var s = initState(m)
    s.livePools[0].add v063LiveAA(1) # compatible PA: prev A, last A
    s.livePools[0].add v063LiveBA(2) # compatible PA: prev B, last A
    s.livePools[0].add v063LiveAB(3) # incompatible PA: last B
    s.monomerN[1] = 10
    s.monomerN0[1] = 10

    check eligiblePoolIndices(m, s, 0) == @[0, 1]
    check livePoolInvariantErrors(m, s).len == 1
    let a = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 2.0 * 10.0
    check abs(a[0] - expected) <= expected * 1.0e-12

    var rng = initRand(6301)
    applyChannel(m, s, 0, rng)

    check s.monomerN[1] == 9
    check s.livePools[2].len == 1
    check poolAcceptsChain(m, 2, s.livePools[2][0])
    check s.livePools[2][0].last == 1
    check s.livePools[2][0].prev == 0
    check s.livePools[0].len == 2
    check eligiblePoolIndices(m, s, 0).len == 1

  test "different-pool termination filters incompatible chains in both pools":
    var m = modelForEngineHardening063()
    var s = initState(m)
    s.livePools[2].add v063LiveAB(10) # compatible PAB
    s.livePools[2].add v063LiveBB(11) # incompatible PAB: prev B, last B
    s.livePools[3].add v063LiveBA(20) # compatible PBA
    s.livePools[3].add v063LiveAA(21) # incompatible PBA: prev A, last A

    let a = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 1.0 * 1.0
    check abs(a[1] - expected) <= expected * 1.0e-12

    var rng = initRand(6302)
    applyChannel(m, s, 1, rng)

    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTermC
    check s.deadChains[0].sequenceText == "A|B|A|B"
    check s.deadChains[0].nMer == @[int32(2), int32(2)]
    check s.livePools[2].len == 1
    check s.livePools[2][0].mers.toText(@["A", "B"]) == "B|B"
    check s.livePools[3].len == 1
    check s.livePools[3][0].mers.toText(@["A", "B"]) == "A|A"
    check livePoolInvariantErrors(m, s).len == 2
