import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc modelForEngineClosure064(): Model =
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
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "PB", kind: pkActive),
    PoolDef(name: "PAB", kind: pkActive),
    PoolDef(name: "PBA", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[0, 1, 1, 0, -1]
  result.poolPenultimateMer = @[-2, -2, 0, 1, -1]
  result.rates = @[
    RateDef(name: "ktr_pa", kConst: 1.0),
    RateDef(name: "ktd_ab_ba", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "transfer_PA_CTA", kind: chMacroTransfer, kId: 0,
               pool1: 0, speciesId: 0, poolOut: 4, speciesOutId: 1),
    KmcChannel(name: "term_d_PAB_PBA", kind: chMacroTermD, kId: 1,
               pool1: 2, pool2: 3, poolOut: 4)
  ]

proc closure064LiveAA(id: int64 = 1): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  result.pushMonomer(0, 100.0)

proc closure064LiveBA(id: int64 = 2): LiveChain =
  result = makeLiveChain(id, "R", 1, 128.0, 2)
  result.pushMonomer(0, 100.0)

proc closure064LiveAB(id: int64 = 3): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  result.pushMonomer(1, 128.0)

proc closure064LiveBB(id: int64 = 4): LiveChain =
  result = makeLiveChain(id, "R", 1, 128.0, 2)
  result.pushMonomer(1, 128.0)

suite "engine phase closure v0.6.4":
  test "transfer from mixed terminal pool filters only by terminal metadata":
    var m = modelForEngineClosure064()
    var s = initState(m)
    s.livePools[0].add closure064LiveAA(1) # compatible PA: prev A, last A
    s.livePools[0].add closure064LiveBA(2) # compatible PA: prev B, last A
    s.livePools[0].add closure064LiveAB(3) # incompatible PA: last B
    s.speciesN[0] = 10

    check eligiblePoolIndices(m, s, 0) == @[0, 1]
    let a = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 2.0 * 10.0
    check abs(a[0] - expected) <= expected * 1.0e-12

    var rng = initRand(6401)
    applyChannel(m, s, 0, rng)

    check s.speciesN[0] == 9
    check s.speciesN[1] == 1
    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTransfer
    check s.deadChains[0].right_end == "CTA"
    check s.deadChains[0].lastMer == 0
    check s.deadChains[0].sequenceText == "A|A" or s.deadChains[0].sequenceText == "B|A"
    check s.livePools[0].len == 2
    check livePoolInvariantErrors(m, s).len == 1

  test "different-pool disproportionation filters incompatible chains in both pools":
    var m = modelForEngineClosure064()
    var s = initState(m)
    s.livePools[2].add closure064LiveAB(10) # compatible PAB
    s.livePools[2].add closure064LiveBB(11) # incompatible PAB: prev B, last B
    s.livePools[3].add closure064LiveBA(20) # compatible PBA
    s.livePools[3].add closure064LiveAA(21) # incompatible PBA: prev A, last A

    let a = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 1.0 * 1.0
    check abs(a[1] - expected) <= expected * 1.0e-12

    var rng = initRand(6402)
    applyChannel(m, s, 1, rng)

    check s.deadChains.len == 2
    var sawAB = false
    var sawBA = false
    var sawH = false
    var sawU = false
    for d in s.deadChains:
      if d.sequenceText == "A|B": sawAB = true
      if d.sequenceText == "B|A": sawBA = true
      if d.formedBy == fbTermD_H and d.right_end == "H": sawH = true
      if d.formedBy == fbTermD_U and d.right_end == "U": sawU = true
    check sawAB
    check sawBA
    check sawH
    check sawU
    check s.livePools[2].len == 1
    check s.livePools[2][0].mers.toText(@["A", "B"]) == "B|B"
    check s.livePools[3].len == 1
    check s.livePools[3][0].mers.toText(@["A", "B"]) == "A|A"
    check livePoolInvariantErrors(m, s).len == 2
