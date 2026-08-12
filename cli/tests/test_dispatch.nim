## Tier 1 tests for slimmc's dispatch logic (dispatch_logic.nim).
##
## These test the sniffing/counting/routing logic in isolation, using
## small temp .model-like files -- no homo/copo engine code is linked
## or invoked here (that end-to-end check is Tier 2, tests/test_dispatch_smoke.py,
## which actually runs the compiled slimmc binary against real models --
## see docs/development/TESTING.md).

import std/[unittest, os, tempfiles]
import ../dispatch_logic

proc writeTempModel(content: string): string =
  let (handle, path) = createTempFile("slimmc_cli_test_", ".model")
  handle.write(content)
  handle.close()
  result = path

suite "monomer counting":
  test "counts exactly one monomer line for a homo-style model":
    let path = writeTempModel("""
desc "one monomer"
monomer M 1.0 100.12
species I 1.0e-5
""")
    defer: removeFile(path)
    check countMonomers(path) == 1

  test "counts two monomer lines for a binary copo-style model":
    let path = writeTempModel("""
desc "binary copolymer"
monomer A 0.2 100.12
monomer B 0.2 128.17
""")
    defer: removeFile(path)
    check countMonomers(path) == 2

  test "counts three monomer lines for a terpolymer copo-style model":
    let path = writeTempModel("""
monomer A 0.2 100.12
monomer B 0.2 128.17
monomer C 0.1 104.15
""")
    defer: removeFile(path)
    check countMonomers(path) == 3

  test "ignores a commented-out monomer line":
    let path = writeTempModel("""
monomer A 0.2 100.12
# monomer B 0.2 128.17
""")
    defer: removeFile(path)
    check countMonomers(path) == 1

  test "ignores monomer keyword appearing only in a comment, not as a directive":
    let path = writeTempModel("""
# this model has one monomer declared below
monomer A 0.2 100.12
""")
    defer: removeFile(path)
    check countMonomers(path) == 1

  test "returns 0 for a model with no monomer declarations":
    let path = writeTempModel("""
desc "no monomers -- a legitimate pure-kinetics model, routes to homo"
species I 1.0e-5
""")
    defer: removeFile(path)
    check countMonomers(path) == 0

  test "returns -1 for a file that does not exist":
    check countMonomers("/this/path/does/not/exist.model") == -1

suite "engine routing":
  test "0 monomers routes to homo (pure-kinetics models have no chain growth at all -- confirmed by homo engine validation cases)":
    check engineForCount(0) == "homo"

  test "1 monomer routes to homo":
    check engineForCount(1) == "homo"

  test "2 monomers routes to copo":
    check engineForCount(2) == "copo"

  test "3 monomers routes to copo":
    check engineForCount(3) == "copo"

  test "4 or more monomers routes to neither (empty string signals dispatcher-level error)":
    check engineForCount(4) == ""
    check engineForCount(10) == ""

echo "[test_dispatch] ok"
