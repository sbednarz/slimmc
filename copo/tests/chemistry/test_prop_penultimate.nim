import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc modelForPenultimatePropTest(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 4
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.species = @[
    SpeciesDef(name: "CTA", c0: 0.0),
    SpeciesDef(name: "Rcta", c0: 0.0)
  ]
  result.pools = @[
    PoolDef(name: "PAB", kind: pkActive),
    PoolDef(name: "PBA", kind: pkActive),
    PoolDef(name: "PBB", kind: pkActive),
    PoolDef(name: "PAA", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[1, 0, 1, 0, -1]
  result.poolPenultimateMer = @[0, 1, 1, 0, -1]
  result.rates = @[
    RateDef(name: "kp_ab_a", kConst: 1.0),
    RateDef(name: "kp_ab_b", kConst: 1.0),
    RateDef(name: "ktc_ab_ba", kConst: 1.0),
    RateDef(name: "ktr_ab", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "prop_PAB_A", kind: chMacroProp, kId: 0,
               pool1: 0, monomerId: 0, poolOut: 1),
    KmcChannel(name: "prop_PAB_B", kind: chMacroProp, kId: 1,
               pool1: 0, monomerId: 1, poolOut: 2),
    KmcChannel(name: "term_c_PAB_PBA", kind: chMacroTermC, kId: 2,
               pool1: 0, pool2: 1, poolOut: 4),
    KmcChannel(name: "transfer_PAB_CTA", kind: chMacroTransfer, kId: 3,
               pool1: 0, speciesId: 0, poolOut: 4, speciesOutId: 1)
  ]

proc liveAB(id: int64 = 1; leftEnd: string = "R"): LiveChain =
  result = makeLiveChain(id, leftEnd, 0, 100.0, 2)
  result.pushMonomer(1, 128.0)

proc liveBA(id: int64 = 2; leftEnd: string = "S"): LiveChain =
  result = makeLiveChain(id, leftEnd, 1, 128.0, 2)
  result.pushMonomer(0, 100.0)

suite "penultimate propagation engine v0.6":
  test "PAB + A moves AB chain to PBA and updates sequence, mass and statistics":
    var m = modelForPenultimatePropTest()
    var s = initState(m)
    s.livePools[0].add liveAB()
    s.monomerN[0] = 1
    s.monomerN0[0] = 1
    var rng = initRand(601)

    applyChannel(m, s, 0, rng)

    check s.monomerN[0] == 0
    check s.livePools[0].len == 0
    check s.livePools[1].len == 1
    let c = s.livePools[1][0]
    check c.dp == 3
    check c.mass == 328.0
    check c.last == 0
    check c.prev == 1
    check c.nMer == @[int32(2), int32(1)]
    check c.mers.toText(@["A", "B"]) == "A|B|A"
    check polymerCompositionCounts(m, s) == @[int64(2), int64(1)]
    check globalDyads(m, s) == @[int64(0), int64(1), int64(1), int64(0)]
    check globalTriads(m, s) == @[int64(0), int64(0), int64(1), int64(0), int64(0), int64(0), int64(0), int64(0)]

  test "PAB + B moves AB chain to PBB":
    var m = modelForPenultimatePropTest()
    var s = initState(m)
    s.livePools[0].add liveAB()
    s.monomerN[1] = 1
    s.monomerN0[1] = 1
    var rng = initRand(602)

    applyChannel(m, s, 1, rng)

    check s.livePools[0].len == 0
    check s.livePools[2].len == 1
    let c = s.livePools[2][0]
    check c.dp == 3
    check c.last == 1
    check c.prev == 1
    check c.nMer == @[int32(1), int32(2)]
    check c.mers.toText(@["A", "B"]) == "A|B|B"

  test "termination by combination works between penultimate pools":
    var m = modelForPenultimatePropTest()
    var s = initState(m)
    s.livePools[0].add liveAB(1, "R")
    s.livePools[1].add liveBA(2, "S")
    var rng = initRand(603)

    applyChannel(m, s, 2, rng)

    check s.livePools[0].len == 0
    check s.livePools[1].len == 0
    check s.deadChains.len == 1
    let d = s.deadChains[0]
    check d.formedBy == fbTermC
    check d.left_end == "R"
    check d.right_end == "S"
    check d.dp == 4
    check d.nMer == @[int32(2), int32(2)]
    check d.sequenceText == "A|B|A|B"

  test "transfer works from penultimate pool without changing monomer inventory":
    var m = modelForPenultimatePropTest()
    var s = initState(m)
    s.livePools[0].add liveAB()
    s.speciesN[0] = 1
    s.speciesN[1] = 0
    s.monomerN = @[int64(2), int64(2)]
    s.monomerN0 = @[int64(3), int64(3)]
    var rng = initRand(604)

    applyChannel(m, s, 3, rng)

    check s.livePools[0].len == 0
    check s.speciesN[0] == 0
    check s.speciesN[1] == 1
    check s.monomerN == @[int64(2), int64(2)]
    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTransfer
    check s.deadChains[0].right_end == "CTA"
    check s.deadChains[0].sequenceText == "A|B"
    check polymerCompositionCounts(m, s) == @[int64(1), int64(1)]
