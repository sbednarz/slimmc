import unittest, random
import copo_types
import copo_kmc
import copo_stats
import copo_sequence

proc baseModel(penultimate: bool): Model =
  result.V = 1.0e-18
  result.dp_max = int64(high(int32))
  result.sequence_mode = "full"
  result.monomers = @[
    MonomerDef(name: "A", c0: 0.0, mw: 100.0),
    MonomerDef(name: "B", c0: 0.0, mw: 128.0)
  ]
  if penultimate:
    result.deadPoolId = 4
    result.pools = @[
      PoolDef(name: "PAA", kind: pkActive), PoolDef(name: "PAB", kind: pkActive),
      PoolDef(name: "PBA", kind: pkActive), PoolDef(name: "PBB", kind: pkActive),
      PoolDef(name: "D", kind: pkDead)
    ]
    result.poolTerminalMer = @[0,1,0,1,-1]
    result.poolPenultimateMer = @[0,0,1,1,-1]
  else:
    result.deadPoolId = 2
    result.pools = @[PoolDef(name:"PA",kind:pkActive), PoolDef(name:"PB",kind:pkActive), PoolDef(name:"D",kind:pkDead)]
    result.poolTerminalMer = @[0,1,-1]
    result.poolPenultimateMer = @[-1,-1,-1]

proc total(xs: seq[int64]): int64 =
  for x in xs: result += x

proc chain2(a,b:int; id:int64=1): LiveChain =
  result = makeLiveChain(id, "R", a, (if a==0:100.0 else:128.0), 2)
  result.pushMonomer(b, (if b==0:100.0 else:128.0))

suite "terminal and penultimate propagation matrix plus microstructure":
  test "all four terminal transitions update terminal, composition and motifs":
    for oldLast in 0..1:
      for added in 0..1:
        var m=baseModel(false)
        m.rates = @[RateDef(name:"kp",kConst:1.0)]
        m.channels = @[KmcChannel(name:"prop",kind:chMacroProp,kId:0,pool1:oldLast,monomerId:added,poolOut:added)]
        var s=initState(m)
        s.livePools[oldLast].add makeLiveChain(1,"R",oldLast,(if oldLast==0:100.0 else:128.0),2)
        s.monomerN[added]=1; s.monomerN0[added]=1
        var rng=initRand(100+oldLast*10+added)
        applyChannel(m,s,0,rng)
        check s.livePools[added].len == 1
        let c=s.livePools[added][0]
        check c.dp == 2
        check c.prev == oldLast
        check c.last == added
        check c.nMer[oldLast] >= 1
        check c.nMer[added] >= 1
        let d=globalDyads(m,s)
        check d[oldLast*2+added] == 1
        check total(d) == 1

  test "all eight penultimate transitions map XY plus Z to YZ and preserve exact motifs":
    for first in 0..1:
      for last in 0..1:
        for added in 0..1:
          var m=baseModel(true)
          let inPool=first*2+last
          let outPool=last*2+added
          m.rates = @[RateDef(name:"kp",kConst:1.0)]
          m.channels = @[KmcChannel(name:"prop",kind:chMacroProp,kId:0,pool1:inPool,monomerId:added,poolOut:outPool)]
          var s=initState(m)
          s.livePools[inPool].add chain2(first,last)
          s.monomerN[added]=1; s.monomerN0[added]=1
          var rng=initRand(300+first*100+last*10+added)
          applyChannel(m,s,0,rng)
          check s.livePools[outPool].len == 1
          let c=s.livePools[outPool][0]
          check c.dp == 3
          check c.prev == last
          check c.last == added
          check c.mers.toText(@["A","B"]) == @["A","B"][first] & "|" & @["A","B"][last] & "|" & @["A","B"][added]
          let dy=globalDyads(m,s)
          check dy[first*2+last] >= 1
          check dy[last*2+added] >= 1
          check total(dy) == 2
          let tr=globalTriads(m,s)
          check tr[first*4+last*2+added] == 1
          check total(tr) == 1

  test "termC reverse join creates the exact boundary dyad and preserves endpoint balance":
    var m=baseModel(true)
    m.rates = @[RateDef(name:"kt",kConst:1.0)]
    m.channels = @[KmcChannel(name:"term",kind:chMacroTermC,kId:0,pool1:1,pool2:2,poolOut:4)]
    var s=initState(m)
    s.livePools[1].add chain2(0,1,1) # AB
    s.livePools[2].add chain2(1,0,2) # BA
    var rng=initRand(999)
    applyChannel(m,s,0,rng)
    check s.deadChains.len == 1
    let d=s.deadChains[0]
    check d.dp == 4
    check d.sequenceText == "A|B|A|B"
    check d.dyads == @[int64(0),int64(2),int64(1),int64(0)]
    check dyadEndpointBalanceOk(d.dyads,2,d.firstMer,d.lastMer)
