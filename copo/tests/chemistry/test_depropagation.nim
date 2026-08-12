import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc modelForDepropTest(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 2
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.pools = @[
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "PB", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[0, 1, -1]
  result.poolPenultimateMer = @[-2, -2, -1]
  result.rates = @[
    RateDef(name: "kdeprop_AA", kConst: 1.0),
    RateDef(name: "kdeprop_BA", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "deprop_PA_to_PA_A", kind: chMacroDeprop, kId: 0,
               pool1: 0, monomerId: 0, poolOut: 0),
    KmcChannel(name: "deprop_PA_to_PB_A", kind: chMacroDeprop, kId: 1,
               pool1: 0, monomerId: 0, poolOut: 1)
  ]

suite "terminal depropagation":
  test "deprop returns terminal monomer and moves shortened chain to previous-terminal pool":
    var m = modelForDepropTest()
    var s = initState(m)
    var c = makeLiveChain(1, "R", 1, 128.0, 2) # B
    c.pushMonomer(0, 100.0)                    # BA, active terminal A -> PA
    s.livePools[0].add c
    s.monomerN[0] = 0
    s.monomerN0[0] = 1
    var rng = initRand(55)

    let a = computePropensities(m, s)
    check a[0] == 0.0
    check a[1] == 1.0
    applyChannel(m, s, 1, rng)

    check s.monomerN[0] == 1
    check s.livePools[0].len == 0
    check s.livePools[1].len == 1
    let shortened = s.livePools[1][0]
    check shortened.dp == 1
    check shortened.mass == 128.0
    check shortened.last == 1
    check shortened.prev == -1
    check shortened.nMer == @[int32(0), int32(1)]
    check shortened.mers.toText(@["A", "B"]) == "B"

  test "deprop is inactive for DP one chains":
    var m = modelForDepropTest()
    var s = initState(m)
    s.livePools[0].add makeLiveChain(1, "R", 0, 100.0, 2)
    let a = computePropensities(m, s)
    check a[0] == 0.0
    check a[1] == 0.0

  test "deprop channel requires matching penultimate terminal":
    var m = modelForDepropTest()
    var s = initState(m)
    var c = makeLiveChain(1, "R", 1, 128.0, 2)
    c.pushMonomer(0, 100.0) # BA, only PA -> PB + A is eligible
    s.livePools[0].add c
    let a = computePropensities(m, s)
    check a[0] == 0.0
    check a[1] == 1.0

  test "two consecutive terminal depropagations update terminal and penultimate":
    var m = modelForDepropTest()
    var s = initState(m)
    var c = makeLiveChain(1, "R", 1, 128.0, 2) # B
    c.pushMonomer(0, 100.0)                    # BA
    c.pushMonomer(0, 100.0)                    # BAA, active terminal A -> PA
    s.livePools[0].add c
    s.monomerN[0] = 0
    var rng = initRand(56)

    let a0 = computePropensities(m, s)
    check a0[0] == 1.0 # ...A-A* -> ...A* + A
    check a0[1] == 0.0
    applyChannel(m, s, 0, rng)

    check s.livePools[0].len == 1
    check s.livePools[1].len == 0
    var afterFirst = s.livePools[0][0]
    check afterFirst.dp == 2
    check afterFirst.last == 0
    check afterFirst.prev == 1
    check afterFirst.mers.toText(@["A", "B"]) == "B|A"
    check s.monomerN[0] == 1

    let a1 = computePropensities(m, s)
    check a1[0] == 0.0
    check a1[1] == 1.0 # ...B-A* -> ...B* + A
    applyChannel(m, s, 1, rng)

    check s.livePools[0].len == 0
    check s.livePools[1].len == 1
    let afterSecond = s.livePools[1][0]
    check afterSecond.dp == 1
    check afterSecond.last == 1
    check afterSecond.prev == -1
    check afterSecond.nMer == @[int32(0), int32(1)]
    check afterSecond.mers.toText(@["A", "B"]) == "B"
    check s.monomerN[0] == 2
    check polymerCompositionCounts(m, s) == @[int64(0), int64(1)]
