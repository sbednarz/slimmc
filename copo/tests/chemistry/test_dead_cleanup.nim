import unittest
import copo_types
import copo_stats

suite "dead summary sequence modes":
  test "composition omits sequence text but preserves statistics":
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    c.pushMonomer(0, 100.0)
    c.pushMonomer(1, 128.0)

    let d = makeDeadSummary(c, "H", fbTermD_H, @["A", "B"], false)
    check d.dp == 4
    check d.nMer == @[int32(2), int32(2)]
    check d.mass == 456.0
    check not d.sequenceStored
    check d.sequenceText == ""
    check d.firstMer == 0
    check d.penultimateMer == 0
    check d.lastMer == 1
    check d.dyads[1] == 2
    check d.dyads[2] == 1

  test "full retains complete sequence text":
    var c = makeLiveChain(1, "R", 0, 100.0, 2)
    c.pushMonomer(1, 128.0)
    let d = makeDeadSummary(c, "H", fbTermD_H, @["A", "B"], true)
    check d.sequenceStored
    check d.sequenceText == "A|B"
