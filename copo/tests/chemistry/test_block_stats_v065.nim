import unittest
import copo_types
import copo_stats
import copo_sequence

proc blockCountOf065(counts: seq[BlockCount]; monomerId, length: int): int64 =
  for bc in counts:
    if bc.monomerId == monomerId and bc.length == int32(length):
      return bc.count
  return 0

proc block065AABBA(): LiveChain =
  result = makeLiveChain(6501, "R", 0, 100.0, 2) # A
  result.pushMonomer(0, 100.0)                    # AA
  result.pushMonomer(1, 128.0)                    # AAB
  result.pushMonomer(1, 128.0)                    # AABB
  result.pushMonomer(0, 100.0)                    # AABBA

suite "exact block statistics v0.6.5":
  test "dead summaries preserve block counts when sequence text is dropped":
    let c = block065AABBA()
    let d = makeDeadSummary(c, "H", fbTermD_H, @["A", "B"], false)

    check d.dp == 5
    check not d.sequenceStored
    check d.sequenceText == ""
    check blockCountOf065(d.blockCounts, 0, 2) == 1 # AA
    check blockCountOf065(d.blockCounts, 0, 1) == 1 # terminal A
    check blockCountOf065(d.blockCounts, 1, 2) == 1 # BB
    check dyadEndpointBalanceOk(d.dyads, 2, d.firstMer, d.lastMer)
    var corruptedDyads = d.dyads
    corruptedDyads[0 * 2 + 1] += 1
    check not dyadEndpointBalanceOk(corruptedDyads, 2, d.firstMer, d.lastMer)

  test "term_c block statistics include boundary merge after reverse join":
    var c1 = makeLiveChain(6502, "R1", 0, 100.0, 2) # A
    c1.pushMonomer(0, 100.0)                         # AA
    c1.pushMonomer(0, 100.0)                         # AAA

    var c2 = makeLiveChain(6503, "R2", 1, 128.0, 2) # B
    c2.pushMonomer(1, 128.0)                         # BB
    c2.pushMonomer(1, 128.0)                         # BBB
    c2.pushMonomer(0, 100.0)                         # BBBA; reverse contributes ABBB

    let d = combineLiveToDead(c1, c2, @["A", "B"], true)

    check d.sequenceText == "A|A|A|A|B|B|B"
    check blockCountOf065(d.blockCounts, 0, 4) == 1
    check blockCountOf065(d.blockCounts, 1, 3) == 1
    check blockCountOf065(d.blockCounts, 0, 3) == 0
    check blockCountOf065(d.blockCounts, 0, 1) == 0

  test "global block counts aggregate live sequences and compact dead summaries":
    var m = Model()
    m.monomers = @[MonomerDef(name: "A", c0: 0.0, mw: 100.0), MonomerDef(name: "B", c0: 0.0, mw: 128.0)]
    m.pools = @[PoolDef(name: "PB", kind: pkActive)]
    var s = State()
    s.livePools = newSeq[seq[LiveChain]](1)

    var live = makeLiveChain(6504, "R", 1, 128.0, 2) # B
    live.pushMonomer(1, 128.0)                       # BB
    s.livePools[0].add live

    var dead = makeDeadSummary(block065AABBA(), "H", fbTermD_H, @["A", "B"], false)
    dead.count = 2
    s.deadChains.add dead

    let counts = globalBlockCounts(m, s)
    check blockCountOf065(counts, 0, 2) == 2 # AA from two compact dead chains
    check blockCountOf065(counts, 0, 1) == 2 # terminal A from two compact dead chains
    check blockCountOf065(counts, 1, 2) == 3 # two dead BB blocks + one live BB block
