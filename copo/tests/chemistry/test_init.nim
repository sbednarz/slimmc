import unittest, random
import copo_types
import copo_kmc
import copo_sequence

proc modelForInitTest(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 2
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.species = @[SpeciesDef(name: "R", c0: 0.0)]
  result.pools = @[
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "PB", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.rates = @[RateDef(name: "ki_a", kConst: 1.0)]
  result.channels = @[
    KmcChannel(name: "init_R_A", kind: chMacroInit, kId: 0,
               speciesId: 0, monomerId: 0, poolOut: 0)
  ]

suite "macro initiation":
  test "init consumes initiator species and monomer, then creates terminal pool chain":
    var m = modelForInitTest()
    var s = initState(m)
    s.speciesN[0] = 1
    s.monomerN[0] = 2
    s.monomerN0[0] = 2
    s.monomerN[1] = 3
    s.monomerN0[1] = 3
    var rng = initRand(11)

    applyChannel(m, s, 0, rng)

    check s.speciesN[0] == 0
    check s.monomerN[0] == 1
    check s.livePools[0].len == 1
    check s.livePools[1].len == 0
    let c = s.livePools[0][0]
    check c.dp == 1
    check c.mass == 100.0
    check c.last == 0
    check c.prev == -1
    check c.nMer == @[int32(1), int32(0)]
    check c.mers.toText(@["A", "B"]) == "A"
