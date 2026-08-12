import unittest, math
import copo_types
import copo_kmc
import copo_stats

proc phaseBModel(): Model =
  result.V = 2.0e-18
  result.dp_max = 100
  result.sequence_mode = "full"
  result.deadPoolId = 2
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  result.species = @[
    SpeciesDef(name: "Q", c0: 0.0),
    SpeciesDef(name: "X", c0: 0.0),
    SpeciesDef(name: "Y", c0: 0.0)
  ]
  result.pools = @[
    PoolDef(name: "PA", kind: pkActive),
    PoolDef(name: "PB", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolTerminalMer = @[0, 1, -1]
  result.poolPenultimateMer = @[-2, -2, -1]
  result.rates = @[
    RateDef(name: "k1", kConst: 2.0),
    RateDef(name: "k2", kConst: 7.0)
  ]

proc chainA(id: int64; dp = 1): LiveChain =
  result = makeLiveChain(id, "R", 0, 100.0, 2)
  for _ in 2 .. dp:
    result.pushMonomer(0, 100.0)

proc closeEnough(actual, expected: float; rtol = 1.0e-12): bool =
  abs(actual - expected) <= max(1.0, abs(expected)) * rtol

suite "copo SSA phase B":
  test "C07 propensity scales with rate, populations, monomer count, and volume":
    var m = phaseBModel()
    m.channels = @[
      KmcChannel(name: "uni", kind: chRxnUni, kId: 0,
                 smallReactants: @[SmallRef(kind: skSpecies, id: 0, stoich: 1)]),
      KmcChannel(name: "bi", kind: chRxnBiDiff, kId: 0,
                 smallReactants: @[SmallRef(kind: skSpecies, id: 0, stoich: 1), SmallRef(kind: skSpecies, id: 1, stoich: 1)]),
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0,
                 pool1: 0, monomerId: 0, poolOut: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0,
                 pool1: 0, monomerId: 0, poolOut: 0)
    ]
    var s = initState(m)
    s.speciesN[0] = 5
    s.speciesN[1] = 3
    s.monomerN[0] = 11
    s.livePools[0] = @[chainA(1, 1), chainA(2, 4), chainA(3, 100)]
    let p = computePropensities(m, s)
    check closeEnough(p[0], 2.0 * 5.0)
    check closeEnough(p[1], 2.0 / (NA * m.V) * 5.0 * 3.0)
    check closeEnough(p[2], 2.0 / (NA * m.V) * 2.0 * 11.0)
    check closeEnough(p[3], 2.0 * 2.0)

    var mRate = m
    mRate.rates[0].kConst *= 3.0
    let pRate = computePropensities(mRate, s)
    for i in 0 ..< p.len:
      check closeEnough(pRate[i], 3.0 * p[i])

    var mVol = m
    mVol.V *= 4.0
    let pVol = computePropensities(mVol, s)
    check closeEnough(pVol[0], p[0])
    check closeEnough(pVol[1], p[1] / 4.0)
    check closeEnough(pVol[2], p[2] / 4.0)
    check closeEnough(pVol[3], p[3])

  test "C07 same-reactant propensity uses n times n-minus-one":
    var m = phaseBModel()
    m.channels = @[
      KmcChannel(name: "same", kind: chRxnBiSame, kId: 1,
                 smallReactants: @[SmallRef(kind: skSpecies, id: 0, stoich: 2)])
    ]
    var s = initState(m)
    s.speciesN[0] = 9
    let actual = computePropensities(m, s)[0]
    let expected = 7.0 / (NA * m.V) * 9.0 * 8.0
    check closeEnough(actual, expected)

  test "C07 unavailable reactants and ineligible chains give zero propensity":
    var m = phaseBModel()
    m.channels = @[
      KmcChannel(name: "uni", kind: chRxnUni, kId: 0,
                 smallReactants: @[SmallRef(kind: skSpecies, id: 0, stoich: 1)]),
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0,
                 pool1: 0, monomerId: 0, poolOut: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0,
                 pool1: 0, monomerId: 0, poolOut: 0)
    ]
    var s = initState(m)
    s.monomerN[0] = 10
    s.livePools[0] = @[chainA(1, 1)]
    let p = computePropensities(m, s)
    check p[0] == 0.0
    check p[1] > 0.0
    check p[2] == 0.0
