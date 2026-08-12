import unittest
import copo_sequence

suite "LinearSequence":
  test "push last prev pop":
    var s = initSequence(0)
    s.push(1)
    s.push(0)
    check s.len == 3
    check s.last == 0
    check s.prev == 1
    check s.pop == 0
    check s.last == 1

  test "appendReverse":
    var a = initSequence(0) # A
    a.push(1)               # AB
    var b = initSequence(1) # B
    b.push(0)               # BA
    a.appendReverse(b)      # ABAB
    check a.data == @[uint8(0), uint8(1), uint8(0), uint8(1)]
