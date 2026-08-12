import unittest, math, tables
import slimmc_types
import slimmc_kmc

proc baseModel(): Model =
  result = initModel()
  result.V = 2.0e-18
  result.T = 298.15
  result.dpMax = 100
  result.species = @[
    SpeciesDef(name: "M", kind: skMonomer, c0: 0.0, mw: 100.0, hasMw: true),
    SpeciesDef(name: "A", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "B", kind: skSpecies, c0: 0.0)
  ]
  for i, sp in result.species:
    result.speciesByName[sp.name] = i
  result.monomerId = 0
  result.pools = @[
    PoolDef(name: "P", kind: pkActive),
    PoolDef(name: "D", kind: pkDead)
  ]
  result.poolByName["P"] = 0
  result.poolByName["D"] = 1
  result.rates = @[
    RateDef(name: "k1", kind: rkFixed, kConst: 2.0),
    RateDef(name: "k2", kind: rkFixed, kConst: 7.0)
  ]
  result.rateByName["k1"] = 0
  result.rateByName["k2"] = 1
  discard result.ensureBuiltinEg("R", 68.0, "test")

proc liveChain(m: Model; dp: int): Chain =
  Chain(eg1: m.egByName["R"], dp: dp, eg2: m.egActive, formedBy: fbInit)

proc stateFor(m: Model; monomer, a, b: int64; live: seq[Chain]): State =
  result = initState(m)
  result.n[0] = monomer
  result.n[1] = a
  result.n[2] = b
  result.pools[0] = live
  result.mExpected = calcMTotal(m, result)
  result.mBalance = result.mExpected

proc closeEnough(actual, expected: float; rtol = 1.0e-12): bool =
  abs(actual - expected) <= max(1.0, abs(expected)) * rtol

suite "homo SSA phase B":
  test "H07 propensity scales exactly with rate, populations, concentration, and volume":
    var m = baseModel()
    m.channels = @[
      KmcChannel(name: "uni", kind: chRxnUni, kId: 0,
                 lhs: @[Term(sp: 1, stoich: 1)]),
      KmcChannel(name: "bi", kind: chRxnBiDiff, kId: 0,
                 lhs: @[Term(sp: 1, stoich: 1), Term(sp: 2, stoich: 1)]),
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0)
    ]
    let s = stateFor(m, 11, 5, 3, @[liveChain(m, 1), liveChain(m, 4), liveChain(m, 100)])
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

  test "H07 same-reactant bimolecular propensity uses n times n-minus-one":
    var m = baseModel()
    m.channels = @[
      KmcChannel(name: "same", kind: chRxnBiSame, kId: 1,
                 lhs: @[Term(sp: 1, stoich: 2)])
    ]
    let s = stateFor(m, 0, 9, 0, @[])
    let p = computePropensities(m, s)[0]
    let expected = 7.0 / (NA * m.V) * 9.0 * 8.0
    check closeEnough(p, expected)

  test "H07 unavailable reactants and ineligible chains give zero propensity":
    var m = baseModel()
    m.channels = @[
      KmcChannel(name: "uni", kind: chRxnUni, kId: 0,
                 lhs: @[Term(sp: 1, stoich: 1)]),
      KmcChannel(name: "prop", kind: chMacroProp, kId: 0, pool1: 0, sp1: 0),
      KmcChannel(name: "deprop", kind: chMacroDeprop, kId: 0, pool1: 0, sp1: 0)
    ]
    let s = stateFor(m, 10, 0, 0, @[liveChain(m, 1)])
    let p = computePropensities(m, s)
    check p[0] == 0.0
    check p[1] > 0.0 # DP=1 can propagate.
    check p[2] == 0.0 # DP=1 cannot depropagate.
