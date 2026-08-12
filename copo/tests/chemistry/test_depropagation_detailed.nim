import unittest, random
import copo_types
import copo_kmc
import copo_stats

proc detailedDepropModel(): Model =
  result.V = 2.0e-18
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
    RateDef(name: "kAA", kConst: 2.0),
    RateDef(name: "kBA", kConst: 3.0),
    RateDef(name: "kAB", kConst: 5.0),
    RateDef(name: "kBB", kConst: 7.0)
  ]
  result.channels = @[
    KmcChannel(name: "AA", kind: chMacroDeprop, kId: 0, pool1: 0, monomerId: 0, poolOut: 0),
    KmcChannel(name: "BA", kind: chMacroDeprop, kId: 1, pool1: 0, monomerId: 0, poolOut: 1),
    KmcChannel(name: "AB", kind: chMacroDeprop, kId: 2, pool1: 1, monomerId: 1, poolOut: 0),
    KmcChannel(name: "BB", kind: chMacroDeprop, kId: 3, pool1: 1, monomerId: 1, poolOut: 1)
  ]

proc seqChain(xs: openArray[int]): LiveChain =
  result = makeLiveChain(1, "R", xs[0], (if xs[0] == 0: 100.0 else: 128.0), 2)
  for i in 1 ..< xs.len:
    result.pushMonomer(xs[i], (if xs[i] == 0: 100.0 else: 128.0))

suite "detailed terminal depropagation":
  test "AA BA AB BB channels select only matching penultimate terminal":
    var m = detailedDepropModel()
    var s = initState(m)
    s.livePools[0].add seqChain([0, 0]) # AA
    s.livePools[0].add seqChain([1, 0]) # BA
    s.livePools[1].add seqChain([0, 1]) # AB
    s.livePools[1].add seqChain([1, 1]) # BB
    let a = computePropensities(m, s)
    check a == @[2.0, 3.0, 5.0, 7.0]

  test "propensity scales with eligible-chain multiplicity and rate":
    var m = detailedDepropModel()
    var s = initState(m)
    for _ in 0 ..< 4: s.livePools[0].add seqChain([0, 0])
    for _ in 0 ..< 3: s.livePools[0].add seqChain([1, 0])
    let a = computePropensities(m, s)
    check a[0] == 8.0
    check a[1] == 9.0
    check a[2] == 0.0
    check a[3] == 0.0

  test "each transition returns its terminal monomer and moves to expected pool":
    for channelId in 0 .. 3:
      var m = detailedDepropModel()
      var s = initState(m)
      let seqs = @[@[0,0], @[1,0], @[0,1], @[1,1]]
      let inputPool = if channelId < 2: 0 else: 1
      let outputPool = if channelId in [0,2]: 0 else: 1
      let returnedMer = if channelId < 2: 0 else: 1
      s.livePools[inputPool].add seqChain(seqs[channelId])
      let before = polymerCompositionCounts(m, s)
      var rng = initRand(100 + channelId)
      applyChannel(m, s, channelId, rng)
      check s.monomerN[returnedMer] == 1
      check s.livePools[inputPool].len + s.livePools[outputPool].len >= 1
      check s.livePools[outputPool].len == 1
      let c = s.livePools[outputPool][0]
      check c.dp == 1
      check c.last == seqs[channelId][0]
      check c.prev == -1
      let after = polymerCompositionCounts(m, s)
      check before[returnedMer] - after[returnedMer] == 1
      check before[1-returnedMer] == after[1-returnedMer]

  test "deprop conserves total A and B units independently":
    var m = detailedDepropModel()
    var s = initState(m)
    s.livePools[0].add seqChain([0, 1, 0])
    s.livePools[1].add seqChain([1, 0, 1])
    let initialPoly = polymerCompositionCounts(m, s)
    let initialFree = s.monomerN
    var rng = initRand(777)
    for _ in 0 ..< 2:
      let a = computePropensities(m, s)
      var chosen = -1
      for i, value in a:
        if value > 0.0: chosen = i; break
      check chosen >= 0
      applyChannel(m, s, chosen, rng)
    let finalPoly = polymerCompositionCounts(m, s)
    for mer in 0 .. 1:
      check initialPoly[mer] + initialFree[mer] == finalPoly[mer] + s.monomerN[mer]

  test "DP one chains contribute zero even with nonzero deprop rates":
    var m = detailedDepropModel()
    var s = initState(m)
    s.livePools[0].add seqChain([0])
    s.livePools[1].add seqChain([1])
    check computePropensities(m, s) == @[0.0, 0.0, 0.0, 0.0]
