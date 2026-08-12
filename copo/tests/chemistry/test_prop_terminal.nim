import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc modelForPropTerminalTest(): Model =
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
  result.rates = @[RateDef(name: "kp_ab", kConst: 1.0)]
  result.channels = @[
    KmcChannel(name: "prop_PA_B", kind: chMacroProp, kId: 0,
               pool1: 0, monomerId: 1, poolOut: 1)
  ]

suite "terminal propagation":
  test "PA + B moves chain to PB and updates sequence, mass and terminal indices":
    var m = modelForPropTerminalTest()
    var s = initState(m)
    s.livePools[0].add makeLiveChain(1, "R", 0, 100.0, 2)
    s.monomerN[1] = 1
    s.monomerN0[1] = 1
    var rng = initRand(22)

    applyChannel(m, s, 0, rng)

    check s.monomerN[1] == 0
    check s.livePools[0].len == 0
    check s.livePools[1].len == 1
    let c = s.livePools[1][0]
    check c.dp == 2
    check c.mass == 228.0
    check c.last == 1
    check c.prev == 0
    check c.nMer == @[int32(1), int32(1)]
    check c.mers.toText(@["A", "B"]) == "A|B"
