import unittest
import copo_sequence

suite "microstructure counts":
  test "AABABBA dyads and triads":
    var s = initSequence(0)
    for x in [0, 1, 0, 1, 1, 0]:
      s.push(x)
    let dy = dyadCounts(s, 2)
    # sequence A A B A B B A
    check dy[0*2+0] == 1 # AA
    check dy[0*2+1] == 2 # AB
    check dy[1*2+0] == 2 # BA
    check dy[1*2+1] == 1 # BB
    let tr = triadCounts(s, 2)
    check tr[(0*2+0)*2+1] == 1 # AAB
    check tr[(0*2+1)*2+0] == 1 # ABA
