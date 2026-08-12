import unittest, random
import copo_types
import copo_kmc
import copo_stats

proc modelForTermDTest(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 1
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.pools = @[
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.rates = @[RateDef(name: "ktd", kConst: 1.0)]
  result.channels = @[
    KmcChannel(name: "term_d_PA_PA", kind: chMacroTermD, kId: 0,
               pool1: 0, pool2: 0, poolOut: 1)
  ]

suite "termination by disproportionation":
  test "termD creates two dead summaries, not one combined chain":
    var m = modelForTermDTest()
    var s = initState(m)
    s.livePools[0].add makeLiveChain(1, "R1", 0, 100.0, 2)
    s.livePools[0].add makeLiveChain(2, "R2", 0, 100.0, 2)
    var rng = initRand(33)

    applyChannel(m, s, 0, rng)

    check s.livePools[0].len == 0
    check s.deadChains.len == 2
    check s.deadChains[0].dp + s.deadChains[1].dp == 2
    let ends = @[s.deadChains[0].right_end, s.deadChains[1].right_end]
    check ends.contains("H")
    check ends.contains("U")
    check s.deadChains[0].formedBy in {fbTermD_H, fbTermD_U}
    check s.deadChains[1].formedBy in {fbTermD_H, fbTermD_U}
