import unittest, os, random
import copo_types
import copo_parser
import copo_kmc
import copo_stats
import copo_sequence

suite "P2 chemistry alignment v0.6.13":
  test "parser accepts elementary rxn, term_x and transfer_m syntax":
    let path = getTempDir() / "slimmc_copo_p2_parser_v0612.model"
    writeFile(path, """
desc "P2 parser"
param kmc_volume 1.6605390671738466e-24
param t_end 0.1
param max_steps 10
monomer A 1.0 100.0
monomer B 1.0 128.0
species PI 1.0
species R 0.0
species Cap 1.0
polymer PA active
polymer PB active
polymer D dead
rate kd 1.0
rate ki 1.0
rate kx 1.0
rate ktrm 1.0
rxn PI -> R + R kd
macro init R + A -> PA ki
macro term_x PA + Cap -> D kx
macro transfer_m PA + A -> D + PA ktrm
""")
    let m = parseModel(path)
    check m.channels[0].kind == chRxnUni
    check m.channels[1].kind == chMacroInit
    check m.channels[2].kind == chMacroTermX
    check m.channels[3].kind == chMacroTransferM

  test "elementary rxn updates species and monomer counters":
    var m: Model
    m.V = 1.0e-18
    m.monomers = @[MonomerDef(name: "A", c0: 0.0, mw: 100.0), MonomerDef(name: "B", c0: 0.0, mw: 128.0)]
    m.species = @[SpeciesDef(name: "PI", c0: 0.0), SpeciesDef(name: "R", c0: 0.0)]
    m.pools = @[PoolDef(name: "D", kind: pkDead)]
    m.poolTerminalMer = @[-1]
    m.poolPenultimateMer = @[-1]
    m.deadPoolId = 0
    m.rates = @[RateDef(name: "kd", kind: rkFixed, kConst: 1.0)]
    m.channels = @[KmcChannel(name: "rxn_PI_to_R_A", kind: chRxnUni, kId: 0,
      smallReactants: @[SmallRef(kind: skSpecies, id: 0, stoich: 1)],
      smallProducts: @[SmallRef(kind: skSpecies, id: 1, stoich: 1),
                       SmallRef(kind: skMonomer, id: 0, stoich: 1)],
      smallDeltas: @[SmallDelta(kind: skSpecies, id: 0, delta: -1),
                     SmallDelta(kind: skSpecies, id: 1, delta: 1),
                     SmallDelta(kind: skMonomer, id: 0, delta: 1)],
      efficiency: 1.0)]
    var s = initState(m)
    s.speciesN = @[int64(1), int64(0)]
    s.monomerN = @[int64(0), int64(0)]
    var rng = initRand(120)

    check computePropensities(m, s)[0] == 1.0
    applyChannel(m, s, 0, rng)

    check s.speciesN == @[int64(0), int64(1)]
    check s.monomerN == @[int64(1), int64(0)]

  test "term_x consumes capping species and creates capped dead summary":
    var m: Model
    m.V = 1.0e-18
    m.dp_max = int64(high(int32))
    m.sequence_mode = "full"
    m.monomers = @[MonomerDef(name: "A", c0: 0.0, mw: 100.0), MonomerDef(name: "B", c0: 0.0, mw: 128.0)]
    m.species = @[SpeciesDef(name: "Cap", c0: 0.0)]
    m.pools = @[PoolDef(name: "PA", kind: pkActive), PoolDef(name: "D", kind: pkDead)]
    m.poolTerminalMer = @[0, -1]
    m.poolPenultimateMer = @[-2, -1]
    m.deadPoolId = 1
    m.rates = @[RateDef(name: "kx", kind: rkFixed, kConst: 1.0)]
    m.channels = @[KmcChannel(name: "term_x_PA_Cap", kind: chMacroTermX, kId: 0,
      pool1: 0, speciesId: 0, poolOut: 1)]
    var s = initState(m)
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0)
    s.livePools[0].add c
    s.speciesN[0] = 1
    var rng = initRand(121)

    applyChannel(m, s, 0, rng)

    check s.speciesN[0] == 0
    check s.livePools[0].len == 0
    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTermX
    check s.deadChains[0].right_end == "Cap"
    check s.deadChains[0].sequenceText == "A|B|A"

  test "transfer_m consumes monomer, terminates chain and creates monomer-born active chain":
    var m: Model
    m.V = 1.0e-18
    m.dp_max = int64(high(int32))
    m.sequence_mode = "full"
    m.monomers = @[MonomerDef(name: "A", c0: 0.0, mw: 100.0), MonomerDef(name: "B", c0: 0.0, mw: 128.0)]
    m.species = newSeq[SpeciesDef]()
    m.pools = @[PoolDef(name: "PA", kind: pkActive), PoolDef(name: "D", kind: pkDead)]
    m.poolTerminalMer = @[0, -1]
    m.poolPenultimateMer = @[-2, -1]
    m.deadPoolId = 1
    m.rates = @[RateDef(name: "ktrm", kind: rkFixed, kConst: 1.0)]
    m.channels = @[KmcChannel(name: "transfer_m_PA_A", kind: chMacroTransferM, kId: 0,
      pool1: 0, monomerId: 0, poolOut: 1, pool2: 0)]
    var s = initState(m)
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0)
    s.livePools[0].add c
    s.monomerN[0] = 1
    s.monomerN0[0] = 1
    var rng = initRand(122)

    applyChannel(m, s, 0, rng)

    check s.monomerN[0] == 0
    check s.deadChains.len == 1
    check s.deadChains[0].formedBy == fbTransferM
    check s.deadChains[0].right_end == "H"
    check s.livePools[0].len == 1
    check s.livePools[0][0].left_end == "A_tr"
    check s.livePools[0][0].dp == 1
    check s.livePools[0][0].last == 0
