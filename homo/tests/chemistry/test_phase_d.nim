import unittest, random, math, tables
import slimmc_types
import slimmc_kmc

proc closeEnough(a, b: float; rtol = 1.0e-12): bool =
  if a == b: return true
  abs(a-b) <= max(abs(a), abs(b))*rtol

proc baseModel(): Model =
  result = initModel()
  result.V = 1.0e-18
  result.T = 298.15
  result.dpMax = 100
  result.species = @[
    SpeciesDef(name: "M", kind: skMonomer, c0: 0.0, mw: 100.0, hasMw: true),
    SpeciesDef(name: "R", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "CTA", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "Rcta", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "CAP", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "A", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "B", kind: skSpecies, c0: 0.0),
    SpeciesDef(name: "C", kind: skSpecies, c0: 0.0)
  ]
  for i, sp in result.species: result.speciesByName[sp.name] = i
  result.monomerId = 0
  result.pools = @[PoolDef(name: "P", kind: pkActive), PoolDef(name: "D", kind: pkDead)]
  result.poolByName["P"] = 0; result.poolByName["D"] = 1
  result.rates = @[RateDef(name: "k", kind: rkFixed, kConst: 2.0)]
  result.rateByName["k"] = 0
  discard result.ensureBuiltinEg("R",68.0,"test")
  discard result.ensureBuiltinEg("Rcta",44.0,"test")
  discard result.ensureBuiltinEg("CAP",40.0,"test")

proc live(m: Model; dp=5): Chain = Chain(eg1: m.egByName["R"], dp: dp, eg2: m.egActive, formedBy: fbInit)
proc state(m: Model; counts: openArray[int64]; living: seq[Chain] = @[]): State =
  result = initState(m)
  for i, n in counts: result.n[i] = n
  result.pools[0] = living
  result.mExpected = calcMTotal(m,result); result.mBalance = result.mExpected

suite "homo chemistry phase D":
  test "H15 macro initiation exact propensity and bookkeeping":
    var m=baseModel()
    m.channels = @[KmcChannel(name: "init",kind: chMacroInit,kId: 0,sp1: 1,sp2: 0,poolOut: 0)]
    var s=state(m,[10'i64,4,0,0,0,0,0,0])
    check closeEnough(computePropensities(m,s)[0],2.0/(NA*m.V)*4.0*10.0)
    let total0=calcMTotal(m,s)
    var rng=initRand(401); applyChannel(m,s,0,rng)
    check s.n[0]==9 and s.n[1]==3
    check s.pools[0].len==1 and s.pools[0][0].dp==1
    check s.pools[0][0].formedBy==fbInit and s.pools[0][0].eg1==m.egByName["R"]
    check calcMTotal(m,s)==total0

  test "H15 initiation requires both radical and monomer":
    var m=baseModel(); m.channels = @[KmcChannel(name: "init",kind: chMacroInit,kId: 0,sp1: 1,sp2: 0,poolOut: 0)]
    check computePropensities(m,state(m,[0'i64,4,0,0,0,0,0,0]))[0]==0
    check computePropensities(m,state(m,[10'i64,0,0,0,0,0,0,0]))[0]==0

  test "H16 transfer-H exact propensity and products":
    var m=baseModel(); m.channels = @[KmcChannel(name: "transfer",kind: chMacroTransferH,kId: 0,pool1: 0,poolOut: 1,sp1: 2,sp2: 3)]
    var s=state(m,[0'i64,0,7,0,0,0,0,0],@[live(m,8),live(m,3)])
    check closeEnough(computePropensities(m,s)[0],2.0/(NA*m.V)*2.0*7.0)
    let total0=calcMTotal(m,s); var rng=initRand(402); applyChannel(m,s,0,rng)
    check s.n[2]==6 and s.n[3]==1
    check s.pools[0].len==1 and s.pools[1].len==1
    check s.pools[1][0].formedBy==fbTransferH and s.pools[1][0].eg2==m.egH
    check calcMTotal(m,s)==total0

  test "H17 term-X exact propensity and cap end group":
    var m=baseModel(); m.channels = @[KmcChannel(name: "termx",kind: chMacroTermX,kId: 0,pool1: 0,poolOut: 1,sp1: 4)]
    var s=state(m,[0'i64,0,0,0,5,0,0,0],@[live(m,6),live(m,9)])
    check closeEnough(computePropensities(m,s)[0],2.0/(NA*m.V)*2.0*5.0)
    let total0=calcMTotal(m,s); var rng=initRand(403); applyChannel(m,s,0,rng)
    check s.n[4]==4 and s.pools[0].len==1 and s.pools[1].len==1
    check s.pools[1][0].formedBy==fbTermX and s.pools[1][0].eg2==m.egByName["CAP"]
    check calcMTotal(m,s)==total0

  test "H18 unimolecular, bimolecular-different, and bimolecular-same propensities":
    var m=baseModel()
    var s=state(m,[0'i64,0,0,0,0,6,4,0])
    m.channels = @[
      KmcChannel(name: "uni",kind: chRxnUni,kId: 0,lhs: @[Term(sp: 5,stoich: 1)],rhs: @[Term(sp: 7,stoich: 1)],eff: 1.0),
      KmcChannel(name: "diff",kind: chRxnBiDiff,kId: 0,lhs: @[Term(sp: 5,stoich: 1),Term(sp: 6,stoich: 1)],rhs: @[Term(sp: 7,stoich: 1)],eff: 1.0),
      KmcChannel(name: "same",kind: chRxnBiSame,kId: 0,lhs: @[Term(sp: 5,stoich: 2)],rhs: @[Term(sp: 7,stoich: 1)],eff: 1.0)
    ]
    let a=computePropensities(m,s)
    check a[0]==12.0
    check closeEnough(a[1],2.0/(NA*m.V)*6.0*4.0)
    check closeEnough(a[2],2.0/(NA*m.V)*6.0*5.0)

  test "H18 elementary reactions apply exact stoichiometry":
    var m=baseModel(); var rng=initRand(404)
    m.channels = @[KmcChannel(name: "same",kind: chRxnBiSame,kId: 0,lhs: @[Term(sp: 5,stoich: 2)],rhs: @[Term(sp: 7,stoich: 1)],eff: 1.0)]
    var s=state(m,[0'i64,0,0,0,0,6,0,0]); applyChannel(m,s,0,rng)
    check s.n[5]==4 and s.n[7]==1
    check s.channelSuccesses[0]==1 and s.channelFailures[0]==0

  test "H19 efficiency zero consumes reactants but forms no products":
    var m=baseModel(); m.channels = @[KmcChannel(name: "eff0",kind: chRxnUni,kId: 0,lhs: @[Term(sp: 5,stoich: 1)],rhs: @[Term(sp: 7,stoich: 1)],eff: 0.0)]
    var s=state(m,[0'i64,0,0,0,0,2,0,0]); var rng=initRand(405); applyChannel(m,s,0,rng)
    check s.n[5]==1 and s.n[7]==0 and s.channelFailures[0]==1

  test "H19 efficiency one always forms products":
    var m=baseModel(); m.channels = @[KmcChannel(name: "eff1",kind: chRxnUni,kId: 0,lhs: @[Term(sp: 5,stoich: 1)],rhs: @[Term(sp: 7,stoich: 1)],eff: 1.0)]
    var s=state(m,[0'i64,0,0,0,0,2,0,0]); var rng=initRand(406); applyChannel(m,s,0,rng)
    check s.n[5]==1 and s.n[7]==1 and s.channelSuccesses[0]==1

  test "H19 intermediate efficiency is statistically correct":
    var m=baseModel(); m.channels = @[KmcChannel(name: "eff",kind: chRxnUni,kId: 0,lhs: @[Term(sp: 5,stoich: 1)],rhs: @[Term(sp: 7,stoich: 1)],eff: 0.25)]
    var successes=0
    const reps=8000
    for seed in 0..<reps:
      var s=state(m,[0'i64,0,0,0,0,1,0,0]); var rng=initRand(10000+seed); applyChannel(m,s,0,rng)
      successes += int(s.channelSuccesses[0])
    check abs(float(successes)/float(reps)-0.25)<0.025
