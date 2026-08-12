import unittest, random, os, tables
import copo_types
import copo_kmc
import copo_stats
import copo_sequence
import copo_parser

proc modelForPoolInvariantTest(): Model =
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
    PoolDef(name: "PAB", kind: pkActive),
    PoolDef(name: "PBA", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[1, 0, -1]
  result.poolPenultimateMer = @[0, 1, -1]
  result.rates = @[
    RateDef(name: "kp", kConst: 1.0),
    RateDef(name: "ktr", kConst: 1.0),
    RateDef(name: "ktc", kConst: 1.0)
  ]
  result.channels = @[
    KmcChannel(name: "prop_PAB_A", kind: chMacroProp, kId: 0,
               pool1: 0, monomerId: 0, poolOut: 1),
    KmcChannel(name: "transfer_PAB_CTA", kind: chMacroTransfer, kId: 1,
               pool1: 0, speciesId: 0, poolOut: 2, speciesOutId: 1),
    KmcChannel(name: "term_c_PAB_PAB", kind: chMacroTermC, kId: 2,
               pool1: 0, pool2: 0, poolOut: 2)
  ]

proc invLiveAB062(id: int64 = 1): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  result.pushMonomer(1, 128.0)

proc invLiveBB062(id: int64 = 2): LiveChain =
  result = makeLiveChain(id, "R", 1, 128.0, 2)
  result.pushMonomer(1, 128.0)

suite "engine pool invariant hardening v0.6.2":
  test "propensity and propagation use invariant-compatible penultimate chains only":
    var m = modelForPoolInvariantTest()
    var s = initState(m)
    s.livePools[0].add invLiveAB062(1)
    s.livePools[0].add invLiveBB062(2) # deliberately wrong for PAB: prev B, last B
    s.monomerN[0] = 10
    s.monomerN0[0] = 10

    let a = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 1.0 * 10.0
    check abs(a[0] - expected) <= expected * 1.0e-12
    check livePoolInvariantErrors(m, s).len == 1

    var rng = initRand(6201)
    applyChannel(m, s, 0, rng)

    check s.livePools[0].len == 1
    check s.livePools[0][0].mers.toText(@["A", "B"]) == "B|B"
    check s.livePools[1].len == 1
    check s.livePools[1][0].mers.toText(@["A", "B"]) == "A|B|A"

  test "transfer from penultimate pool ignores chains with stale terminal/penultimate metadata":
    var m = modelForPoolInvariantTest()
    var s = initState(m)
    s.livePools[0].add invLiveBB062(2) # incompatible with PAB
    s.speciesN[0] = 10
    let a0 = computePropensities(m, s)
    check a0[1] == 0.0

    s.livePools[0].add invLiveAB062(1)
    let a1 = computePropensities(m, s)
    let expected = 1.0 / (NA * m.V) * 1.0 * 10.0
    check abs(a1[1] - expected) <= expected * 1.0e-12

  test "same-pool termination counts only compatible chains":
    var m = modelForPoolInvariantTest()
    var s = initState(m)
    s.livePools[0].add invLiveAB062(1)
    s.livePools[0].add invLiveBB062(2)
    check computePropensities(m, s)[2] == 0.0
    s.livePools[0].add invLiveAB062(3)
    let expected = 1.0 / (NA * m.V) * 2.0 * 1.0
    check abs(computePropensities(m, s)[2] - expected) <= expected * 1.0e-12

  test "parser infers penultimate pool metadata independent of macro order":
    let path = getTempDir() / "slimmc_copo_parser_v062_order.model"
    writeFile(path, """
desc "parser v0.6.2 order-independent pool metadata"
param kmc_volume 1.0e-19
param t_end 1.0
param max_steps 10
param seed 1
monomer A 0.1 100.0
monomer B 0.1 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer PAB active
polymer PBA active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_ab_a 30.0
rate kp_pa_b 20.0
macro prop PAB + A -> PBA kp_ab_a
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + B -> PAB kp_pa_b
""")
    let m = parseModel(path)
    let pab = m.poolByName["PAB"]
    let pba = m.poolByName["PBA"]
    check m.poolTerminalMer[pab] == 1
    check m.poolPenultimateMer[pab] == 0
    check m.poolTerminalMer[pba] == 0
    check m.poolPenultimateMer[pba] == 1
