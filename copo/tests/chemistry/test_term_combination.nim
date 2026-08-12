import unittest
import copo_types
import copo_stats

suite "termination by combination":
  test "termC stores seq1 followed by reverse(seq2) and preserves composition":
    var c1 = makeLiveChain(1, "R1", 0, 100.0, 2) # A
    c1.pushMonomer(1, 128.0)                       # AB
    var c2 = makeLiveChain(2, "R2", 1, 128.0, 2) # B
    c2.pushMonomer(0, 100.0)                       # BA

    let d = combineLiveToDead(c1, c2, @["A", "B"], true)

    check d.formedBy == fbTermC
    check d.dp == 4
    check d.mass == 456.0
    check d.left_end == "R1"
    check d.right_end == "R2"
    check d.nMer == @[int32(2), int32(2)]
    check d.sequenceStored
    check d.sequenceText == "A|B|A|B"
    check d.firstMer == 0
    check d.penultimateMer == 0
    check d.lastMer == 1
    check d.dyads[0 * 2 + 1] == 2 # AB
    check d.dyads[1 * 2 + 0] == 1 # BA
    check dyadEndpointBalanceOk(d.dyads, 2, d.firstMer, d.lastMer)
