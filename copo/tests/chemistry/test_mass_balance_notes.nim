import unittest
import copo_types
import copo_stats

suite "chemistry invariants v0.2":
  test "monomer inventory can be reconstructed from free monomer plus live/dead polymer counts":
    var m = Model()
    m.monomers = @[
      MonomerDef(name: "A", c0: 0.0, mw: 100.0),
      MonomerDef(name: "B", c0: 0.0, mw: 128.0)
    ]
    m.pools = @[PoolDef(name: "PA", kind: pkActive), PoolDef(name: "D", kind: pkDead)]
    m.deadPoolId = 1
    var s = State()
    s.monomerN0 = @[int64(3), int64(3)]
    s.monomerN = @[int64(1), int64(2)]
    s.livePools = newSeq[seq[LiveChain]](2)
    var live = makeLiveChain(1, "R", 0, 100.0, 2)
    live.pushMonomer(1, 128.0)
    s.livePools[0].add live
    s.deadChains.add makeDeadSummary(makeLiveChain(2, "R", 0, 100.0, 2), "H", fbTermD_H, @["A", "B"], true)

    let pc = polymerCompositionCounts(m, s)
    check pc == @[int64(2), int64(1)]
    check s.monomerN[0] + pc[0] == s.monomerN0[0]
    check s.monomerN[1] + pc[1] == s.monomerN0[1]

  test "live pool name matches terminal mer after terminal propagation convention":
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    check c.last == 1
    check c.prev == 0
