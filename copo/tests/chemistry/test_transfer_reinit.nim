import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc modelForTransferTest(): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.deadPoolId = 2
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
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[0, 1, -1]
  result.poolPenultimateMer = @[-2, -2, -1]
  result.rates = @[
    RateDef(name: "ktr_a", kConst: 1.0),
    RateDef(name: "ki_cta_b", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "transfer_PA_CTA", kind: chMacroTransfer, kId: 0,
               pool1: 0, speciesId: 0, poolOut: 2, speciesOutId: 1),
    KmcChannel(name: "init_Rcta_B", kind: chMacroInit, kId: 1,
               speciesId: 1, monomerId: 1, poolOut: 1)
  ]

suite "transfer and reinitiation":
  test "transfer creates dead summary and a new radical species":
    var m = modelForTransferTest()
    var s = initState(m)
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0) # ABA, terminal A -> PA
    s.livePools[0].add c
    s.speciesN[0] = 1
    s.speciesN[1] = 0
    var rng = initRand(44)

    applyChannel(m, s, 0, rng)

    check s.speciesN[0] == 0
    check s.speciesN[1] == 1
    check s.livePools[0].len == 0
    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTransfer
    check s.deadChains[0].right_end == "CTA"
    check s.deadChains[0].dp == 3
    check s.deadChains[0].nMer == @[int32(2), int32(1)]
    check s.deadChains[0].sequenceText == "A|B|A"

  test "radical made by transfer reinits through ordinary macro init":
    var m = modelForTransferTest()
    var s = initState(m)
    s.speciesN[1] = 1
    s.monomerN[1] = 2
    s.monomerN0[1] = 2
    var rng = initRand(45)

    applyChannel(m, s, 1, rng)

    check s.speciesN[1] == 0
    check s.monomerN[1] == 1
    check s.livePools[1].len == 1
    let c = s.livePools[1][0]
    check c.left_end == "Rcta"
    check c.dp == 1
    check c.last == 1
    check c.mers.toText(@["A", "B"]) == "B"

  test "transfer preserves mass composition microstructure even when sequence text is dropped":
    var m = modelForTransferTest()
    m.dp_max = int64(high(int32))
    m.sequence_mode = "composition"
    var s = initState(m)
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0) # ABA, terminal A -> PA
    s.livePools[0].add c
    s.speciesN[0] = 1
    var rng = initRand(46)

    applyChannel(m, s, 0, rng)

    check s.deadChains.len == 1
    let d = s.deadChains[0]
    check d.dp == 3
    check d.mass == 328.0
    check d.nMer == @[int32(2), int32(1)]
    check d.firstMer == 0
    check d.lastMer == 0
    check d.sequenceStored == false
    check d.sequenceText == ""
    check d.dyads == @[int64(0), int64(1), int64(1), int64(0)]
    check d.triads == @[int64(0), int64(0), int64(1), int64(0), int64(0), int64(0), int64(0), int64(0)]
    check polymerCompositionCounts(m, s) == @[int64(2), int64(1)]

  test "transfer does not change monomer inventory, reinitiation consumes exactly one monomer":
    var m = modelForTransferTest()
    var s = initState(m)
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0) # ABA, terminal A -> PA
    s.livePools[0].add c
    s.monomerN0 = @[int64(3), int64(3)]
    s.monomerN = @[int64(2), int64(2)]
    s.speciesN[0] = 1
    var rng = initRand(47)

    applyChannel(m, s, 0, rng)
    check s.monomerN == @[int64(2), int64(2)]
    check polymerCompositionCounts(m, s) == @[int64(2), int64(1)]

    applyChannel(m, s, 1, rng)
    check s.monomerN == @[int64(2), int64(1)]
    check polymerCompositionCounts(m, s) == @[int64(2), int64(2)]
